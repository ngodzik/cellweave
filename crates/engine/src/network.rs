//! A protein network compiled from a specification, and run with RK4.
//!
//! [`cellweave_core::NetworkSpec`] is the hypothesis as written; [`Network`] is
//! it made runnable. Compiling resolves every name once, so the hot loop walks
//! indices and never a string, and a name that resolves to nothing is an error
//! at load rather than a silent zero in the middle of a run.

use std::sync::Arc;

use cellweave_core::network::{INPUT_NAMES, Term};
use cellweave_core::{
    CellInputs, CellMechanics, CellweaveError, Concentration, MechanicsSpec, NetworkSpec,
    ProteinState, SignalingNetwork, TimeStep, Volume,
};

/// Where a term reads its source from, once resolved.
#[derive(Debug, Clone, Copy)]
enum Source {
    Oxygen,
    Egf,
    Crowding,
    Stress,
    Node(usize),
}

/// The source an environment input reads, by the name a file gives it.
/// `None` for a name that is not an input.
fn input_source(name: &str) -> Option<Source> {
    match name {
        "oxygen" => Some(Source::Oxygen),
        "egf" => Some(Source::Egf),
        "crowding" => Some(Source::Crowding),
        "stress" => Some(Source::Stress),
        _ => None,
    }
}

/// A term with its source resolved and its shape reduced to what the loop needs.
#[derive(Debug, Clone, Copy)]
enum Compiled {
    /// `rate * f(source) * (1 - self)`, `f` being Hill when `hill` is set.
    ActivatedBy {
        source: Source,
        rate: f64,
        hill: Option<(f64, f64)>,
    },
    /// `rate * (1 - source) * (1 - self)`.
    ActivatedByAbsenceOf { source: Source, rate: f64 },
    /// `rate * source`.
    ProducedBy { source: Source, rate: f64 },
    /// `- rate * self`.
    Decays { rate: f64 },
    /// `- rate * (self - value)`.
    RelaxesTo { value: f64, rate: f64 },
    /// `- rate * max(0, (threshold - source) / threshold)`.
    DamagedBelow {
        source: Source,
        threshold: f64,
        rate: f64,
    },
}

impl Compiled {
    /// The rate this term moves its node at, which is what a mutation changes.
    fn rate_mut(&mut self) -> &mut f64 {
        match self {
            Compiled::ActivatedBy { rate, .. }
            | Compiled::ActivatedByAbsenceOf { rate, .. }
            | Compiled::ProducedBy { rate, .. }
            | Compiled::Decays { rate }
            | Compiled::RelaxesTo { rate, .. }
            | Compiled::DamagedBelow { rate, .. } => rate,
        }
    }

    /// The kind of term, as a network file names it.
    fn kind(&self) -> &'static str {
        match self {
            Compiled::ActivatedBy { .. } => "activated_by",
            Compiled::ActivatedByAbsenceOf { .. } => "activated_by_absence_of",
            Compiled::ProducedBy { .. } => "produced_by",
            Compiled::Decays { .. } => "decays",
            Compiled::RelaxesTo { .. } => "relaxes_to",
            Compiled::DamagedBelow { .. } => "damaged_below",
        }
    }
}

/// How far one mutation can move a rate: up to this factor either way.
///
/// Large on purpose. A mutation that changes a rate by a percent is invisible
/// against the noise of the lattice, and selection can only act on what it can
/// see; doubling or halving one rate is the scale of a gene switched up or
/// down.
pub const MUTATION_SPREAD: f64 = 2.0;

/// The rate a gained link starts at, before the same spread a scaled rate
/// gets: a weak link, of the order of the network's own terms, that has to be
/// strengthened by later mutations to matter much.
pub const GAINED_RATE: f64 = 0.1;

/// One heritable change to a network.
#[derive(Debug, Clone, PartialEq)]
pub struct Mutation {
    /// The node whose terms changed, by name.
    pub node: String,
    /// Which of its terms, counting from zero in the order the file lists
    /// them; for a gained link, the position it was added at, the last.
    pub term: usize,
    /// What kind of term that is, as the file names it.
    pub kind: &'static str,
    /// What happened to it.
    pub change: Change,
}

/// What a mutation did to a term.
#[derive(Debug, Clone, PartialEq)]
pub enum Change {
    /// Its rate was multiplied by this factor, between `1 / MUTATION_SPREAD`
    /// and [`MUTATION_SPREAD`].
    Scaled(f64),
    /// It is new: the node now also reads `source`, at `rate`.
    Gained {
        /// The node or input the new term reads.
        source: String,
        /// The rate it starts at.
        rate: f64,
    },
    /// It was removed.
    Lost,
}

impl std::fmt::Display for Mutation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match &self.change {
            Change::Scaled(factor) => {
                write!(
                    f,
                    "{} {} #{} x{:.2}",
                    self.node, self.kind, self.term, factor
                )
            }
            Change::Gained { source, rate } => {
                write!(
                    f,
                    "{} gains {} {} @{:.2}",
                    self.node, self.kind, source, rate
                )
            }
            Change::Lost => write!(f, "{} loses {} #{}", self.node, self.kind, self.term),
        }
    }
}

/// A protein network, ready to run.
///
/// Cheap to clone, which is what a dividing cell does: a daughter is a copy of
/// its parent's state and of the same compiled terms.
#[derive(Debug, Clone)]
pub struct Network {
    /// Current value of each node, in the spec's order.
    state: ProteinState,
    /// Name of each node, same order, shared by every copy.
    names: Arc<[String]>,
    /// Terms per node, same order.
    terms: Vec<Vec<Compiled>>,
    /// Index of the growth, survival and hypoxia nodes, and of the
    /// differentiation and apoptosis nodes when the network has them.
    growth: usize,
    survival: usize,
    hypoxia: usize,
    differentiation: Option<usize>,
    apoptosis: Option<usize>,
    /// How the growth node becomes mechanics.
    mechanics: MechanicsSpec,
    /// Volume this lineage divides against.
    reference_volume: u32,
    /// Scratch for RK4, so a step allocates nothing: the state as plain floats,
    /// the midpoint, and the four stages.
    x: Vec<f64>,
    mid: Vec<f64>,
    k: [Vec<f64>; 4],
}

impl Network {
    /// Compile a specification into something runnable.
    ///
    /// # Errors
    ///
    /// [`CellweaveError::Config`] when the network has no nodes, when two nodes
    /// share a name, when a term or a role names something that does not exist,
    /// or when a Hill exponent or half-maximal level is not positive.
    pub fn from_spec(spec: &NetworkSpec, reference_volume: u32) -> Result<Self, CellweaveError> {
        if spec.nodes.is_empty() {
            return Err(CellweaveError::Config("a network needs a node".into()));
        }

        let roles = spec.roles.as_ref().ok_or_else(|| {
            CellweaveError::Config(
                "the network declares no roles, and inherits none: say which node is growth, \
                 which is survival, and which is the hypoxia response"
                    .into(),
            )
        })?;
        let mechanics = spec.mechanics.ok_or_else(|| {
            CellweaveError::Config(
                "the network declares no mechanics, and inherits none: say how the growth node \
                 becomes a target volume"
                    .into(),
            )
        })?;

        let mut index = std::collections::HashMap::new();
        for (i, node) in spec.nodes.iter().enumerate() {
            if INPUT_NAMES.contains(&node.name.as_str()) {
                return Err(CellweaveError::Config(format!(
                    "node `{}` takes the name of an environment input; inputs are {:?}",
                    node.name, INPUT_NAMES
                )));
            }
            if index.insert(node.name.clone(), i).is_some() {
                return Err(CellweaveError::Config(format!(
                    "two nodes are named `{}`",
                    node.name
                )));
            }
        }

        let resolve = |name: &str, whose: &str| -> Result<Source, CellweaveError> {
            if let Some(input) = input_source(name) {
                return Ok(input);
            }
            index.get(name).copied().map(Source::Node).ok_or_else(|| {
                let mut known: Vec<&str> = INPUT_NAMES.to_vec();
                known.extend(spec.nodes.iter().map(|n| n.name.as_str()));
                CellweaveError::Config(format!(
                    "node `{whose}` reads `{name}`, which is neither a node nor an input; known names are {known:?}"
                ))
            })
        };

        let mut terms = Vec::with_capacity(spec.nodes.len());
        for node in &spec.nodes {
            let mut compiled = Vec::with_capacity(node.terms.len());
            for term in &node.terms {
                compiled.push(match term {
                    Term::ActivatedBy { source, rate, hill } => {
                        let hill = match hill {
                            Some(h) if h.n <= 0.0 || h.k <= 0.0 => {
                                return Err(CellweaveError::Config(format!(
                                    "node `{}` has a Hill term with n = {} and k = {}; both must be positive",
                                    node.name, h.n, h.k
                                )));
                            }
                            Some(h) => Some((h.n, h.k)),
                            None => None,
                        };
                        Compiled::ActivatedBy {
                            source: resolve(source, &node.name)?,
                            rate: *rate,
                            hill,
                        }
                    }
                    Term::ActivatedByAbsenceOf { source, rate } => {
                        Compiled::ActivatedByAbsenceOf {
                            source: resolve(source, &node.name)?,
                            rate: *rate,
                        }
                    }
                    Term::ProducedBy { source, rate } => Compiled::ProducedBy {
                        source: resolve(source, &node.name)?,
                        rate: *rate,
                    },
                    Term::Decays { rate } => Compiled::Decays { rate: *rate },
                    Term::RelaxesTo { value, rate } => Compiled::RelaxesTo {
                        value: *value,
                        rate: *rate,
                    },
                    Term::DamagedBelow {
                        source,
                        threshold,
                        rate,
                    } => {
                        if *threshold <= 0.0 {
                            return Err(CellweaveError::Config(format!(
                                "node `{}` is damaged below {threshold}, which must be positive",
                                node.name
                            )));
                        }
                        Compiled::DamagedBelow {
                            source: resolve(source, &node.name)?,
                            threshold: *threshold,
                            rate: *rate,
                        }
                    }
                });
            }
            terms.push(compiled);
        }

        let role = |name: &str, which: &str| -> Result<usize, CellweaveError> {
            index.get(name).copied().ok_or_else(|| {
                CellweaveError::Config(format!(
                    "the {which} role names `{name}`, which is not a node"
                ))
            })
        };

        let mut state = ProteinState::zeros(spec.nodes.len());
        for (i, node) in spec.nodes.iter().enumerate() {
            state.values[i] = Concentration(node.initial.clamp(0.0, 1.0));
        }

        Ok(Self {
            state,
            names: spec.nodes.iter().map(|n| n.name.clone()).collect(),
            terms,
            growth: role(&roles.growth, "growth")?,
            survival: role(&roles.survival, "survival")?,
            hypoxia: role(&roles.hypoxia_response, "hypoxia response")?,
            differentiation: roles
                .differentiation
                .as_deref()
                .map(|name| role(name, "differentiation"))
                .transpose()?,
            apoptosis: roles
                .apoptosis
                .as_deref()
                .map(|name| role(name, "apoptosis"))
                .transpose()?,
            mechanics,
            reference_volume,
            x: vec![0.0; spec.nodes.len()],
            mid: vec![0.0; spec.nodes.len()],
            k: std::array::from_fn(|_| vec![0.0; spec.nodes.len()]),
        })
    }

    /// Reference volume of this lineage, which the target volume and the
    /// division threshold are both measured against.
    ///
    /// A daughter inherits it unchanged, so the size at which cells divide is a
    /// property of the lineage rather than of each cell's own birth size, which
    /// is how a cell type keeps a characteristic size across generations.
    pub fn reference_volume(&self) -> u32 {
        self.reference_volume
    }

    /// The names of the nodes, in the order the specification lists them,
    /// which is also the order of a snapshot's protein levels.
    pub fn node_names(&self) -> &[String] {
        &self.names
    }

    /// The node the specification named as survival. A cell whose survival
    /// collapses is necrotic.
    pub fn survival(&self) -> f64 {
        self.state.values[self.survival].0
    }

    /// The node the specification named as the hypoxia response. A cell whose
    /// response is high has stopped cycling: hypoxia arrests the cell cycle
    /// well before it kills.
    pub fn hypoxia_response(&self) -> f64 {
        self.state.values[self.hypoxia].0
    }

    /// Scale one rate, chosen by `which`, by a factor set by `how`, and say what
    /// changed. Both are uniform draws in `[0, 1)`, taken by the caller so that
    /// the network holds no random generator of its own.
    ///
    /// Every term of every node is equally likely, and the factor is uniform
    /// in its logarithm, so a mutation is as likely to halve a rate as to
    /// double it. Only the rate changes, never the shape of a term: a mutation
    /// turns a gene up or down, it does not rewire the network.
    ///
    /// `None` for a network without a single term, which has nothing a
    /// mutation could change.
    pub fn mutate(&mut self, which: f64, how: f64) -> Option<Mutation> {
        let count: usize = self.terms.iter().map(Vec::len).sum();
        let mut pick = ((which * count as f64) as usize).min(count.saturating_sub(1));
        let factor = MUTATION_SPREAD.powf(2.0 * how - 1.0);
        for (node, terms) in self.terms.iter_mut().enumerate() {
            if pick < terms.len() {
                let term = &mut terms[pick];
                *term.rate_mut() *= factor;
                return Some(Mutation {
                    node: self.names[node].clone(),
                    term: pick,
                    kind: term.kind(),
                    change: Change::Scaled(factor),
                });
            }
            pick -= terms.len();
        }
        None
    }

    /// Change the network's wiring rather than a rate, and say what changed.
    /// `draws` are five uniform draws in `[0, 1)`, taken by the caller.
    ///
    /// Half the time a term is lost: one of all the network's terms, every
    /// one equally likely, is removed. The other half a link is gained: a
    /// node, every one equally likely, starts to read a source, any input of
    /// the environment or any node, itself included, either switched on by it
    /// (`activated_by`) or by its absence (`activated_by_absence_of`), at
    /// [`GAINED_RATE`] spread as a scaled rate is. That is how a gene comes
    /// to answer a signal it ignored, or stops answering one.
    ///
    /// `None` when there is nothing to lose, a network without a single term.
    pub fn rewire(&mut self, draws: [f64; 5]) -> Option<Mutation> {
        let [mode, which, source, sign, how] = draws;
        let pick = |n: usize, u: f64| ((u * n as f64) as usize).min(n.saturating_sub(1));
        if mode < 0.5 {
            let node = pick(self.terms.len(), which);
            let sources = INPUT_NAMES.len() + self.names.len();
            let s = pick(sources, source);
            let (source, source_name) = match INPUT_NAMES.get(s) {
                Some(&name) => (input_source(name)?, name.to_owned()),
                None => {
                    let n = s - INPUT_NAMES.len();
                    (Source::Node(n), self.names[n].clone())
                }
            };
            let rate = GAINED_RATE * MUTATION_SPREAD.powf(2.0 * how - 1.0);
            let term = if sign < 0.5 {
                Compiled::ActivatedBy {
                    source,
                    rate,
                    hill: None,
                }
            } else {
                Compiled::ActivatedByAbsenceOf { source, rate }
            };
            let kind = term.kind();
            self.terms[node].push(term);
            return Some(Mutation {
                node: self.names[node].clone(),
                term: self.terms[node].len() - 1,
                kind,
                change: Change::Gained {
                    source: source_name,
                    rate,
                },
            });
        }
        let count: usize = self.terms.iter().map(Vec::len).sum();
        let mut at = pick(count, which);
        for (node, terms) in self.terms.iter_mut().enumerate() {
            if at < terms.len() {
                let lost = terms.remove(at);
                return Some(Mutation {
                    node: self.names[node].clone(),
                    term: at,
                    kind: lost.kind(),
                    change: Change::Lost,
                });
            }
            at -= terms.len();
        }
        None
    }

    /// Make a recorded mutation again: scale the same rate by the same
    /// factor, lose the same term, or gain the same link at the same rate.
    /// Replaying every mutation of a lineage in order, on the network it
    /// started from, gives back the lineage's network exactly.
    ///
    /// `None`, with nothing changed, when the mutation does not fit: a node
    /// this network does not have, or a term that is not where the mutation
    /// says, which happens when an earlier mutation it relied on, a lost or a
    /// gained term on the same node, was left out.
    pub(crate) fn apply(&mut self, mutation: &Mutation) -> Option<()> {
        let node = self.names.iter().position(|n| *n == mutation.node)?;
        match &mutation.change {
            Change::Scaled(factor) => {
                let term = self.terms[node]
                    .get_mut(mutation.term)
                    .filter(|t| t.kind() == mutation.kind)?;
                *term.rate_mut() *= factor;
            }
            Change::Lost => {
                self.terms[node]
                    .get(mutation.term)
                    .filter(|t| t.kind() == mutation.kind)?;
                self.terms[node].remove(mutation.term);
            }
            Change::Gained { source, rate } => {
                let source = match input_source(source) {
                    Some(input) => input,
                    None => Source::Node(self.names.iter().position(|n| n == source)?),
                };
                self.terms[node].push(if mutation.kind == "activated_by" {
                    Compiled::ActivatedBy {
                        source,
                        rate: *rate,
                        hill: None,
                    }
                } else {
                    Compiled::ActivatedByAbsenceOf {
                        source,
                        rate: *rate,
                    }
                });
            }
        }
        Some(())
    }

    /// Take over another network's protein levels, when it has as many
    /// nodes: a cell whose wiring changes keeps the state it is in.
    pub(crate) fn take_state(&mut self, from: &Network) {
        if from.state.values.len() == self.state.values.len() {
            self.state = from.state.clone();
        }
    }

    /// The node the specification named as differentiation, or `None` when it
    /// named none. Its level is the chance per step that the cell commits to
    /// differentiate.
    pub fn differentiation(&self) -> Option<f64> {
        self.differentiation.map(|i| self.state.values[i].0)
    }

    /// The node the specification named as apoptosis, or `None` when it named
    /// none. Its level is the chance per step that the cell starts its own
    /// programmed death.
    pub fn apoptosis(&self) -> Option<f64> {
        self.apoptosis.map(|i| self.state.values[i].0)
    }

    /// The chance per step this network would give of starting apoptosis,
    /// were its cell to sit in `inputs` for `time` units of network time from
    /// where it stands now. Worked out on a copy: the network itself does not
    /// move. `None` for a network without an apoptosis role.
    ///
    /// How a cell would respond somewhere it is not is what separates a
    /// mutation that changes a response from one that only changes where the
    /// cell happens to be.
    pub fn apoptosis_under(&self, inputs: &CellInputs, time: u32) -> Option<f64> {
        self.apoptosis?;
        let mut probe = self.clone();
        for _ in 0..time {
            probe.step(TimeStep(1.0), inputs);
        }
        probe.apoptosis()
    }

    /// The value of one source, given the current state and the environment.
    fn read(state: &[f64], inputs: &CellInputs, source: Source) -> f64 {
        match source {
            Source::Oxygen => inputs.o2.0.clamp(0.0, 1.0),
            Source::Egf => inputs.egf.0,
            Source::Crowding => inputs.crowding.0.clamp(0.0, 1.0),
            Source::Stress => inputs.stress.0.max(0.0),
            Source::Node(i) => state[i],
        }
    }

    /// Rate of change of every node, into `out`.
    ///
    /// Free of `self` so that the caller can hold the terms and a scratch buffer
    /// at once, which is what keeps a step from allocating.
    fn derivatives(terms: &[Vec<Compiled>], state: &[f64], inputs: &CellInputs, out: &mut [f64]) {
        for (i, node_terms) in terms.iter().enumerate() {
            let own = state[i];
            let mut rate = 0.0;
            for term in node_terms {
                rate += match *term {
                    Compiled::ActivatedBy {
                        source,
                        rate: k,
                        hill,
                    } => {
                        let s = Self::read(state, inputs, source);
                        let driver = match hill {
                            None => s,
                            Some((n, half)) => {
                                let sn = s.max(0.0).powf(n);
                                sn / (half.powf(n) + sn)
                            }
                        };
                        k * driver * (1.0 - own)
                    }
                    Compiled::ActivatedByAbsenceOf { source, rate: k } => {
                        k * (1.0 - Self::read(state, inputs, source)) * (1.0 - own)
                    }
                    Compiled::ProducedBy { source, rate: k } => {
                        k * Self::read(state, inputs, source)
                    }
                    Compiled::Decays { rate: k } => -k * own,
                    Compiled::RelaxesTo { value, rate: k } => -k * (own - value),
                    Compiled::DamagedBelow {
                        source,
                        threshold,
                        rate: k,
                    } => {
                        let s = Self::read(state, inputs, source);
                        -k * ((threshold - s) / threshold).max(0.0)
                    }
                };
            }
            out[i] = rate;
        }
    }

    /// Classic RK4 over `dt`.
    ///
    /// Four evaluations, the middle two at the half step, blended one, two, two,
    /// one. Only the result is clamped to the unit range, as the intermediate
    /// stages are slopes rather than states.
    fn rk4(&mut self, dt: f64, inputs: &CellInputs) {
        let n = self.state.values.len();
        for i in 0..n {
            self.x[i] = self.state.values[i].0;
        }

        Self::derivatives(&self.terms, &self.x, inputs, &mut self.k[0]);
        for stage in 0..3 {
            let step = if stage < 2 { 0.5 * dt } else { dt };
            for i in 0..n {
                self.mid[i] = self.x[i] + step * self.k[stage][i];
            }
            let (done, todo) = self.k.split_at_mut(stage + 1);
            debug_assert_eq!(done.len(), stage + 1);
            Self::derivatives(&self.terms, &self.mid, inputs, &mut todo[0]);
        }

        for i in 0..n {
            let delta = (dt / 6.0)
                * (self.k[0][i] + 2.0 * self.k[1][i] + 2.0 * self.k[2][i] + self.k[3][i]);
            self.state.values[i] = Concentration((self.x[i] + delta).clamp(0.0, 1.0));
        }
    }
}

impl SignalingNetwork for Network {
    fn step(&mut self, dt: TimeStep, inputs: &CellInputs) {
        self.rk4(dt.0, inputs);
    }

    fn state(&self) -> &ProteinState {
        &self.state
    }

    fn mechanics(&self) -> CellMechanics {
        let growth = self.state.values[self.growth].0;
        let m = self.mechanics;
        let ratio = 1.0 + m.growth_coefficient * growth;
        let ratio = m.max_target.map_or(ratio, |max| ratio.min(max));
        let target = (f64::from(self.reference_volume) * ratio).max(0.0) as u32;
        CellMechanics {
            target_volume: Volume(target),
            lambda_volume: m.lambda_volume,
            j_medium: m.j_medium,
            j_self: m.j_self_base + m.j_self_from_growth * (1.0 - growth),
        }
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used, clippy::panic)]

    use super::*;
    use cellweave_core::network::{Hill, NodeSpec, Roles};

    /// The network that used to be written as literals in the source.
    const RAS_ERK: &str = include_str!("../../../examples/networks/ras-erk.toml");

    fn ras_erk() -> Network {
        let spec: NetworkSpec = toml::from_str(RAS_ERK).expect("the built-in network parses");
        Network::from_spec(&spec, 200).expect("the built-in network compiles")
    }

    fn at(o2: f64, egf: f64) -> CellInputs {
        CellInputs {
            o2: Concentration(o2),
            egf: Concentration(egf),
            vegf: Concentration(0.0),
            crowding: Concentration(0.0),
            stress: Concentration(0.0),
        }
    }

    fn after(steps: u32, o2: f64, egf: f64) -> Vec<f64> {
        let mut net = ras_erk();
        for _ in 0..steps {
            net.step(TimeStep(1.0), &at(o2, egf));
        }
        net.state().values.iter().map(|c| c.0).collect()
    }

    /// Numbers captured from the hardcoded network before it was replaced.
    ///
    /// Bit for bit, not approximately: the terms are summed in the order the
    /// file lists them, which is the order the old expressions added them, so
    /// there is no reason for a single digit to differ. A change here means the
    /// configurable network is not the same hypothesis, which is exactly what
    /// this is here to catch.
    #[test]
    fn the_built_in_network_matches_the_one_it_replaced() {
        #[rustfmt::skip]
        let cases: [(u32, f64, f64, [f64; 5]); 9] = [
            (1,  0.8,  0.2, [0.115920000000, 0.087273010004, 0.116228133667, 0.090589121667, 0.692551847772]),
            (5,  0.8,  0.2, [0.222916797074, 0.263169229210, 0.187040004493, 0.306151690709, 0.773072487888]),
            (40, 0.8,  0.2, [0.249999977385, 0.333333273023, 0.193548387094, 0.451608501700, 1.000000000000]),
            (1,  0.05, 0.2, [0.115920000000, 0.087273010004, 0.362365362104, 0.188937568542, 0.552470990480]),
            (5,  0.05, 0.2, [0.222916797074, 0.263169229210, 0.530068481696, 0.881237056294, 0.294661809665]),
            (40, 0.05, 0.2, [0.249999977385, 0.333333273023, 0.532710280374, 1.000000000000, 0.188725322011]),
            (1,  1.0,  0.0, [0.037041875000, 0.060105748792, 0.030338541667, 0.060503541667, 0.686250342282]),
            (5,  1.0,  0.0, [0.011157959788, 0.040001480605, 0.004112382360, 0.035817480783, 0.637188154950]),
            (40, 1.0,  0.0, [0.000000307531, 0.000002417035, 0.000000000105, 0.000001383521, 0.501136552845]),
        ];

        for (steps, o2, egf, expected) in cases {
            let got = after(steps, o2, egf);
            for (i, want) in expected.iter().enumerate() {
                assert!(
                    (got[i] - want).abs() < 1e-12,
                    "after {steps} steps at oxygen {o2} and EGF {egf}, node {i} is {:.12} and was {want:.12}",
                    got[i]
                );
            }
        }
    }

    #[test]
    fn the_mechanics_match_the_ones_they_replaced() {
        let mut net = ras_erk();
        for _ in 0..40 {
            net.step(TimeStep(1.0), &at(0.8, 0.2));
        }
        let m = net.mechanics();
        assert_eq!(m.target_volume, Volume(299));
        assert_eq!(m.lambda_volume, 50.0);
        assert_eq!(m.j_medium, 16.0);
        assert!(
            (m.j_self - 7.333333815816).abs() < 1e-11,
            "j_self is {}",
            m.j_self
        );
    }

    fn one_node(terms: Vec<Term>) -> NetworkSpec {
        NetworkSpec {
            nodes: vec![NodeSpec {
                name: "a".into(),
                initial: 0.0,
                terms,
            }],
            extends: None,
            roles: Some(Roles {
                growth: "a".into(),
                survival: "a".into(),
                hypoxia_response: "a".into(),
                differentiation: None,
                apoptosis: None,
            }),
            mechanics: Some(MechanicsSpec {
                growth_coefficient: 1.0,
                lambda_volume: 1.0,
                j_medium: 1.0,
                j_self_base: 1.0,
                j_self_from_growth: 0.0,
                max_target: None,
            }),
        }
    }

    #[test]
    fn replaying_a_lineage_gives_back_its_network() {
        let mut lineage = ras_erk();
        let mut replayed = ras_erk();
        let draws = [
            [0.1, 0.3, 0.6, 0.2, 0.9],
            [0.7, 0.5, 0.1, 0.8, 0.4],
            [0.2, 0.9, 0.95, 0.6, 0.1],
            [0.6, 0.05, 0.3, 0.3, 0.7],
        ];
        let mut made = Vec::new();
        for d in draws {
            made.extend(lineage.rewire(d));
            made.extend(lineage.mutate(d[1], d[4]));
        }
        assert!(made.len() >= 6, "{made:?}");
        for m in &made {
            assert!(replayed.apply(m).is_some(), "{m}");
        }
        assert_eq!(format!("{lineage:?}"), format!("{replayed:?}"));
    }

    #[test]
    fn a_ceiling_caps_the_target_however_high_the_growth() {
        let mut spec: NetworkSpec = toml::from_str(RAS_ERK).unwrap();
        if let Some(m) = spec.mechanics.as_mut() {
            m.max_target = Some(1.2);
        }
        let mut capped = Network::from_spec(&spec, 200).unwrap();
        let mut free = ras_erk();
        for _ in 0..40 {
            capped.step(TimeStep(1.0), &at(0.8, 0.2));
            free.step(TimeStep(1.0), &at(0.8, 0.2));
        }
        assert_eq!(free.mechanics().target_volume, Volume(299));
        assert_eq!(capped.mechanics().target_volume, Volume(240));
    }

    #[test]
    fn a_term_reading_a_name_that_is_neither_node_nor_input_is_refused() {
        let spec = one_node(vec![Term::ActivatedBy {
            source: "p53".into(),
            rate: 1.0,
            hill: None,
        }]);

        let refused = Network::from_spec(&spec, 100);

        let Err(CellweaveError::Config(message)) = refused else {
            panic!("an unknown source must be refused");
        };
        assert!(
            message.contains("p53"),
            "the message must name it: {message}"
        );
        assert!(
            message.contains("oxygen"),
            "and say what is known: {message}"
        );
    }

    #[test]
    fn a_role_naming_something_that_is_not_a_node_is_refused() {
        let mut spec = one_node(vec![]);
        spec.roles.as_mut().unwrap().survival = "bcl2".into();

        let refused = Network::from_spec(&spec, 100);

        let Err(CellweaveError::Config(message)) = refused else {
            panic!("an unknown role must be refused");
        };
        assert!(
            message.contains("survival") && message.contains("bcl2"),
            "{message}"
        );
    }

    #[test]
    fn two_nodes_with_one_name_are_refused() {
        let mut spec = one_node(vec![]);
        spec.nodes.push(spec.nodes[0].clone());

        assert!(matches!(
            Network::from_spec(&spec, 100),
            Err(CellweaveError::Config(_))
        ));
    }

    #[test]
    fn a_node_that_takes_an_input_name_is_refused() {
        let mut spec = one_node(vec![]);
        spec.nodes[0].name = "oxygen".into();
        spec.roles = Some(Roles {
            growth: "oxygen".into(),
            survival: "oxygen".into(),
            hypoxia_response: "oxygen".into(),
            differentiation: None,
            apoptosis: None,
        });

        let refused = Network::from_spec(&spec, 100);

        let Err(CellweaveError::Config(message)) = refused else {
            panic!("shadowing an input must be refused");
        };
        assert!(message.contains("oxygen"), "{message}");
    }

    #[test]
    fn an_empty_network_is_refused() {
        let mut spec = one_node(vec![]);
        spec.nodes.clear();
        assert!(Network::from_spec(&spec, 100).is_err());
    }

    #[test]
    fn a_hill_term_with_an_impossible_shape_is_refused() {
        for hill in [Hill { n: 0.0, k: 0.5 }, Hill { n: 2.0, k: 0.0 }] {
            let spec = one_node(vec![Term::ActivatedBy {
                source: "egf".into(),
                rate: 1.0,
                hill: Some(hill),
            }]);
            assert!(
                Network::from_spec(&spec, 100).is_err(),
                "n = {} and k = {} must be refused",
                hill.n,
                hill.k
            );
        }
    }

    #[test]
    fn each_term_does_what_its_name_says() {
        // Read the slope, not the integrated value: over one very short step the
        // change divided by the step is the derivative, which is exactly what
        // each term's documentation states. The step has to be short because the
        // second order error of the scheme is about the curvature times half the
        // step, and the curvature here is of order one.
        let dt = 1e-6;
        let cases: [(&str, Term, f64, f64, f64); 6] = [
            // rate * source * (1 - self), with source 0.5 and self 0.25
            (
                "activated_by",
                Term::ActivatedBy {
                    source: "egf".into(),
                    rate: 2.0,
                    hill: None,
                },
                0.25,
                0.5,
                2.0 * 0.5 * 0.75,
            ),
            // rate * (1 - source) * (1 - self)
            (
                "activated_by_absence_of",
                Term::ActivatedByAbsenceOf {
                    source: "oxygen".into(),
                    rate: 1.0,
                },
                0.2,
                0.25,
                1.0 * 0.75 * 0.8,
            ),
            // rate * source, no saturation: self plays no part
            (
                "produced_by",
                Term::ProducedBy {
                    source: "egf".into(),
                    rate: 2.0,
                },
                0.9,
                0.5,
                2.0 * 0.5,
            ),
            // - rate * self
            ("decays", Term::Decays { rate: 0.5 }, 0.8, 0.0, -0.5 * 0.8),
            // - rate * (self - value), here pulling up from below the set point
            (
                "relaxes_to",
                Term::RelaxesTo {
                    value: 0.4,
                    rate: 2.0,
                },
                0.1,
                0.0,
                -2.0 * (0.1 - 0.4),
            ),
            // - rate * (threshold - source) / threshold, at a third of the way down
            (
                "damaged_below",
                Term::DamagedBelow {
                    source: "oxygen".into(),
                    threshold: 0.3,
                    rate: 3.0,
                },
                0.5,
                0.2,
                -3.0 * (0.3 - 0.2) / 0.3,
            ),
        ];

        for (name, term, initial, input, expected_slope) in cases {
            let source_is_oxygen = matches!(term.source(), Some("oxygen"));
            let mut spec = one_node(vec![term]);
            spec.nodes[0].initial = initial;
            let mut net = Network::from_spec(&spec, 100).unwrap();
            let inputs = if source_is_oxygen {
                at(input, 0.0)
            } else {
                at(0.0, input)
            };

            net.step(TimeStep(dt), &inputs);

            let slope = (net.state().values[0].0 - initial) / dt;
            assert!(
                (slope - expected_slope).abs() < 1e-5,
                "{name}: slope {slope:.9}, documented {expected_slope:.9}"
            );
        }
    }

    #[test]
    fn a_node_at_zero_cannot_be_damaged_below_zero() {
        let mut spec = one_node(vec![Term::DamagedBelow {
            source: "oxygen".into(),
            threshold: 0.5,
            rate: 10.0,
        }]);
        spec.nodes[0].initial = 0.0;
        let mut net = Network::from_spec(&spec, 100).unwrap();

        for _ in 0..50 {
            net.step(TimeStep(1.0), &at(0.0, 0.0));
        }

        assert_eq!(net.state().values[0].0, 0.0);
    }

    #[test]
    fn decay_pulls_a_node_down_and_a_hill_term_switches_more_sharply() {
        // Decay on its own, from a node that starts high.
        let mut spec = one_node(vec![Term::Decays { rate: 1.0 }]);
        spec.nodes[0].initial = 1.0;
        let mut net = Network::from_spec(&spec, 100).unwrap();
        net.step(TimeStep(0.1), &at(0.0, 0.0));
        assert!(net.state().values[0].0 < 1.0 && net.state().values[0].0 > 0.85);

        // A Hill term with a high exponent is nearly off below its half point
        // and nearly on above it, where a proportional term is neither.
        let hill = |n: f64| {
            let spec = one_node(vec![Term::ActivatedBy {
                source: "egf".into(),
                rate: 10.0,
                hill: Some(Hill { n, k: 0.5 }),
            }]);
            let run = |egf: f64| {
                let mut net = Network::from_spec(&spec, 100).unwrap();
                for _ in 0..200 {
                    net.step(TimeStep(0.1), &at(0.0, egf));
                }
                net.state().values[0].0
            };
            (run(0.3), run(0.7))
        };
        let (soft_low, soft_high) = hill(1.0);
        let (sharp_low, sharp_high) = hill(8.0);
        assert!(
            sharp_low < soft_low,
            "a sharp term must be quieter below the half point"
        );
        assert!(sharp_high > soft_high, "and louder above it");
    }

    #[test]
    fn a_daughter_is_a_copy_that_then_lives_its_own_life() {
        let mut parent = ras_erk();
        for _ in 0..10 {
            parent.step(TimeStep(1.0), &at(0.8, 0.2));
        }
        let mut daughter = parent.clone();
        assert_eq!(parent.state().values, daughter.state().values);

        // Same inputs, same trajectory; different inputs, different fate.
        for _ in 0..30 {
            parent.step(TimeStep(1.0), &at(0.8, 0.2));
            daughter.step(TimeStep(1.0), &at(0.02, 0.2));
        }
        assert!(parent.survival() > 0.9);
        assert!(daughter.survival() < 0.25);
        assert_eq!(parent.reference_volume(), daughter.reference_volume());
    }

    #[test]
    fn a_mutation_scales_one_rate_and_the_parent_keeps_its_own() {
        let parent = ras_erk();
        let mut daughter = parent.clone();
        // The last of the eleven terms is BCL2's damage below the anoxia
        // threshold; a draw near one doubles it.
        let m = daughter.mutate(0.99, 0.999_999).unwrap();
        assert_eq!(
            (m.node.as_str(), m.term, m.kind),
            ("bcl2", 3, "damaged_below")
        );
        let Change::Scaled(factor) = m.change else {
            panic!("a rate was scaled: {m}")
        };
        assert!((factor - MUTATION_SPREAD).abs() < 1e-5);

        let survive = |mut net: Network| {
            for _ in 0..40 {
                net.step(TimeStep(1.0), &at(0.05, 0.2));
            }
            net.survival()
        };
        assert!(
            survive(daughter) < survive(parent),
            "doubled damage kills faster in the dark, and only in the daughter"
        );
    }

    #[test]
    fn a_draw_in_the_middle_changes_nothing() {
        let mut net = ras_erk();
        let m = net.mutate(0.0, 0.5).unwrap();
        assert_eq!(
            (m.node.as_str(), m.term, m.kind),
            ("ras", 0, "activated_by")
        );
        assert_eq!(m.change, Change::Scaled(1.0));
        for _ in 0..5 {
            net.step(TimeStep(1.0), &at(0.8, 0.2));
        }
        let mut reference = ras_erk();
        for _ in 0..5 {
            reference.step(TimeStep(1.0), &at(0.8, 0.2));
        }
        assert_eq!(net.state().values, reference.state().values);
    }

    #[test]
    fn a_crowded_cell_reads_its_crowding_and_answers_with_its_differentiation_node() {
        let mut spec = one_node(vec![
            Term::ActivatedBy {
                source: "crowding".into(),
                rate: 1.0,
                hill: None,
            },
            Term::Decays { rate: 1.0 },
        ]);
        if let Some(roles) = spec.roles.as_mut() {
            roles.differentiation = Some("a".into());
        }
        let mut net = Network::from_spec(&spec, 100).unwrap();
        let squeezed = CellInputs {
            crowding: Concentration(0.6),
            ..at(1.0, 0.0)
        };
        for _ in 0..40 {
            net.step(TimeStep(1.0), &squeezed);
        }
        // Steady state of k c (1 - a) = a: a = c / (1 + c).
        let level = net.differentiation().unwrap();
        assert!((level - 0.6 / 1.6).abs() < 1e-6, "{level}");

        if let Some(roles) = spec.roles.as_mut() {
            roles.differentiation = Some("nowhere".into());
        }
        assert!(
            Network::from_spec(&spec, 100).is_err(),
            "a role naming no node"
        );
        assert_eq!(
            ras_erk().differentiation(),
            None,
            "a tumour network has none"
        );
    }

    #[test]
    fn a_gained_link_makes_a_node_answer_what_it_ignored() {
        let mut net = ras_erk();
        // Gain (below one half), on the first node, RAS; from the first
        // source, oxygen; switched on by it; at the middle of the spread.
        let m = net.rewire([0.0, 0.0, 0.0, 0.0, 0.5]).unwrap();
        assert_eq!(m.to_string(), "ras gains activated_by oxygen @0.10");
        assert_eq!(m.term, 2, "added after RAS's two terms");
        // With no EGF at all, the healthy RAS stays off; this one rises with
        // the oxygen it now reads.
        let ras = |mut n: Network| {
            for _ in 0..40 {
                n.step(TimeStep(1.0), &at(1.0, 0.0));
            }
            n.state().values[0].0
        };
        assert!(ras(ras_erk()) < 0.01);
        assert!(ras(net) > 0.1);
    }

    #[test]
    fn a_lost_term_takes_its_answer_away() {
        let mut net = ras_erk();
        // Loss (one half or above), of the very first term: RAS switched on
        // by EGF.
        let m = net.rewire([0.9, 0.0, 0.0, 0.0, 0.0]).unwrap();
        assert_eq!(m.to_string(), "ras loses activated_by #0");
        assert_eq!(m.change, Change::Lost);
        let ras = |mut n: Network| {
            for _ in 0..40 {
                n.step(TimeStep(1.0), &at(0.8, 0.2));
            }
            n.state().values[0].0
        };
        assert!(ras(ras_erk()) > 0.2, "EGF switches the healthy RAS on");
        assert!(ras(net) < 0.01, "and nothing switches this one on");
    }
}
