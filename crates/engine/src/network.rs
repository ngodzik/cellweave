//! A protein network compiled from a specification, and run with RK4.
//!
//! [`cellweave_core::NetworkSpec`] is the hypothesis as written; [`Network`] is
//! it made runnable. Compiling resolves every name once, so the hot loop walks
//! indices and never a string, and a name that resolves to nothing is an error
//! at load rather than a silent zero in the middle of a run.

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
    Node(usize),
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

/// A protein network, ready to run.
///
/// Cheap to clone, which is what a dividing cell does: a daughter is a copy of
/// its parent's state and of the same compiled terms.
#[derive(Debug, Clone)]
pub struct Network {
    /// Current value of each node, in the spec's order.
    state: ProteinState,
    /// Terms per node, same order.
    terms: Vec<Vec<Compiled>>,
    /// Index of the growth, survival and hypoxia nodes.
    growth: usize,
    survival: usize,
    hypoxia: usize,
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
            match name {
                "oxygen" => Ok(Source::Oxygen),
                "egf" => Ok(Source::Egf),
                other => index.get(other).copied().map(Source::Node).ok_or_else(|| {
                    let mut known: Vec<&str> = INPUT_NAMES.to_vec();
                    known.extend(spec.nodes.iter().map(|n| n.name.as_str()));
                    CellweaveError::Config(format!(
                        "node `{whose}` reads `{other}`, which is neither a node nor an input; known names are {known:?}"
                    ))
                }),
            }
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
            terms,
            growth: role(&spec.roles.growth, "growth")?,
            survival: role(&spec.roles.survival, "survival")?,
            hypoxia: role(&spec.roles.hypoxia_response, "hypoxia response")?,
            mechanics: spec.mechanics,
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

    /// The value of one source, given the current state and the environment.
    fn read(state: &[f64], inputs: &CellInputs, source: Source) -> f64 {
        match source {
            Source::Oxygen => inputs.o2.0.clamp(0.0, 1.0),
            Source::Egf => inputs.egf.0,
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
        let target = (f64::from(self.reference_volume) * (1.0 + m.growth_coefficient * growth))
            .max(0.0) as u32;
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
            drugs: vec![],
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
            roles: Roles {
                growth: "a".into(),
                survival: "a".into(),
                hypoxia_response: "a".into(),
            },
            mechanics: MechanicsSpec {
                growth_coefficient: 1.0,
                lambda_volume: 1.0,
                j_medium: 1.0,
                j_self_base: 1.0,
                j_self_from_growth: 0.0,
            },
        }
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
        spec.roles.survival = "bcl2".into();

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
        spec.roles = Roles {
            growth: "oxygen".into(),
            survival: "oxygen".into(),
            hypoxia_response: "oxygen".into(),
        };

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
}
