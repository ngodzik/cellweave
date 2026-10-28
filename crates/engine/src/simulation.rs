//! The loop that steps the three scales together.
//!
//! It lives in the engine rather than in the binary so that whatever runs a
//! simulation, the command line or a viewer, runs the same one. Two copies of
//! the loop would be two models that drift apart without a test noticing.
//!
//! The lifecycle is a typestate: a [`Simulation<Configured>`] has its tissue
//! seeded but its field not yet solved, and only [`Simulation::start`], which
//! solves it, hands back a [`Simulation<Running>`] that can be stepped. Stepping
//! against an unsettled field is a compile error rather than a quiet artefact.

use std::marker::PhantomData;

use cellweave_core::{
    CellId, CellInputs, CellKind, CellSnapshot, CellweaveError, CloneId, Concentration, Exposure,
    NetworkSpec, OxygenSource, Pos2, Seeding, SignalingNetwork, SimSnapshot, TimeStep, Volume,
};
use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};

use crate::coupling::{
    apply_background_uptake, apply_uptake, hold_medium_bath, hold_obstacles_at, mean_per_cell,
};
use crate::cpm::{CpmLattice, OBSTACLE};
use crate::evolution::CloneRecord;
use crate::network::Network;
use crate::pde::ScalarField;

/// The seed every run uses unless told otherwise.
///
/// Fixed, so that two runs of one configuration are the same run byte for
/// byte, which is what lets a difference between two configurations be
/// attributed to what was changed.
pub const DEFAULT_SEED: u64 = 42;

/// A cell divides once it reaches this multiple of its lineage's reference volume.
const DIVISION_RATIO: f64 = 1.4;

/// Monte Carlo steps between two checks for cells ready to divide, in a tissue
/// grown from one cell.
///
/// Dividing scans the whole lattice, so it is not something to attempt on every
/// step, and a cell needs time to grow back anyway. A sheet checks every step
/// instead: there a cell reaches the threshold only when a neighbour leaves,
/// and checking in batches would make the whole layer divide in waves, every
/// cell freed in twenty steps splitting at once.
const DIVISION_CHECK_EVERY: u64 = 20;

/// Refuse to go past this many living cells, as a guard against runaway
/// growth. Counted among the living rather than among every cell ever born,
/// since a tissue that turns over keeps its size while its births add up.
const MAX_CELLS: usize = 2000;

/// Oxygen diffusion coefficient, in pixels squared per unit time.
const OXYGEN_DIFFUSION: f64 = 0.2;

/// The field is relaxed until no pixel moves by more than this fraction of the
/// uptake per step.
///
/// Relative on purpose: the per step change is what diffusion has not yet
/// balanced against the uptake, so a tolerance above the uptake itself would
/// declare the field settled without having balanced anything.
const OXYGEN_TOLERANCE: f64 = 0.01;

/// Give up on relaxing the field after this many steps, loudly.
const OXYGEN_MAX_STEPS: usize = 20_000;

/// Before the first step the field is solved from scratch, and a per-step
/// tolerance is a poor guard against a slow front creeping out from a small
/// source: each step moves less than the tolerance while the total drifts for
/// hundreds of steps. Once, at the start, the field is settled far more
/// tightly so that every later step is an increment on a true equilibrium.
const OXYGEN_WARMUP_TOLERANCE: f64 = 1e-8;

/// And it is allowed to take its time doing so.
const OXYGEN_WARMUP_MAX_STEPS: usize = 2_000_000;

/// Units of network time that pass in one lattice step, taken in sub steps of
/// one unit each so the RK4 scheme stays well inside its stability limit.
///
/// The three clocks have to agree. The field is relaxed to equilibrium every
/// lattice step, so it is instantaneous. A cell grows to division in about
/// twenty lattice steps. The network settles in about fifteen units of its own
/// time, so at one unit per lattice step it would answer barely faster than the
/// cell cycle, and a cell's fate would trail its surroundings by most of a
/// generation. Eight units per step brings the answer down to two lattice
/// steps: field, then network, then cycle, each well ahead of the next, which is
/// the order biology has them in even if the ratios here are compressed.
const NETWORK_TIME_PER_MCS: u32 = 8;

/// Survival signal below which a cell is necrotic.
const DEATH_THRESHOLD: f64 = 0.25;

/// Hypoxia response above which a cell stops dividing while staying alive.
const QUIESCENCE_THRESHOLD: f64 = 0.45;

/// Pixels a necrotic cell gives up per lattice step until nothing is left.
///
/// Dead tissue is cleared, slowly, and its place taken by whatever grows into
/// it. Without this a tumour is bounded by its box; with it, by turnover.
const RESORPTION_RATE: u32 = 2;

/// Pixels a differentiated cell gives up per lattice step as it leaves the
/// layer.
///
/// Faster than resorption: leaving is something a living cell does, moving up
/// out of the basal layer, not debris waiting to be cleared. A cell of a
/// hundred pixels is gone in about seventeen steps, under a generation.
const SHEDDING_RATE: u32 = 6;

/// Pixels an apoptotic cell gives up per lattice step as it breaks up.
///
/// Faster still than shedding: apoptosis packs a cell into fragments that
/// its neighbours clear within hours, where necrotic debris lingers for days.
/// A cell of two hundred pixels is gone in about twenty-five steps.
const FRAGMENTATION_RATE: u32 = 8;

/// The lethal stress diffuses like oxygen. How far it reaches into the tissue
/// is set by how much the cells take up, which the exposure states, so one
/// diffusion coefficient for both keeps the comparison between the two fields
/// a comparison of their uptakes alone.
const STRESS_DIFFUSION: f64 = OXYGEN_DIFFUSION;

/// What a simulation is asked to be: the box, its temperature, where the oxygen
/// comes from, and the seed.
///
/// The hypothesis inside each cell is not here; it is handed to
/// [`Simulation::new`] separately, because it is the one thing an experiment
/// varies while holding the rest fixed.
#[derive(Debug, Clone, PartialEq)]
pub struct Scenario {
    /// Grid width in pixels.
    pub width: u32,
    /// Grid height in pixels.
    pub height: u32,
    /// CPM temperature.
    pub temperature: f64,
    /// Where the oxygen comes from, which also decides where the first cell goes.
    pub oxygen_source: OxygenSource,
    /// What one consuming pixel takes up per unit time, as a fraction of the
    /// source level.
    pub oxygen_uptake: f64,
    /// Seed of the only random generator in the run.
    pub seed: u64,
    /// What the box holds at the start: one tumour cell, or a sheet of normal
    /// tissue.
    pub seeding: Seeding,
    /// The lethal stress the tissue is exposed to, if any. `None` leaves the
    /// stress field unsolved and every cell reading zero, so a run without
    /// exposure is the run it was before exposure existed.
    pub exposure: Option<Exposure>,
    /// Clones whose founding mutation is, in this run, drawn but not made.
    ///
    /// Everything else happens as it would have: the same numbers are drawn,
    /// the clone is founded with the same identifier and the same parent, and
    /// its cells are the cells they would have been. Only their network is
    /// their parent's, unchanged. Until that difference makes a difference,
    /// the run is the run without the suppression, bit for bit; a mutation
    /// that never makes one leaves it so to the end.
    pub suppress: Vec<CloneId>,
    /// Probability that a daughter is born with one rate of its network
    /// changed, founding a clone of its own. Zero keeps every daughter an
    /// exact copy, and draws nothing from the generator for it, so a run
    /// without mutation is the same run it was before mutation existed.
    pub mutation_rate: f64,
}

/// Lifecycle stage of a [`Simulation`]. Sealed: the stages are the ones below.
pub trait SimState: sealed::Sealed {}

/// Tissue seeded, field not yet solved.
#[derive(Debug)]
pub struct Configured;

/// Field solved; the simulation can be stepped.
#[derive(Debug)]
pub struct Running;

impl SimState for Configured {}
impl SimState for Running {}

mod sealed {
    pub trait Sealed {}
    impl Sealed for super::Configured {}
    impl Sealed for super::Running {}
}

/// What a cell is doing, as the loop decides it.
///
/// The one place the classification lives, so that what divides in the loop
/// and what anything outside it shows as dividing cannot disagree.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Fate {
    /// Alive, and divides once it has grown large enough.
    Cycling,
    /// Alive and holding its place, but not dividing: hypoxia has arrested its
    /// cell cycle.
    Quiescent,
    /// Dead, and giving up its pixels until nothing is left.
    Necrotic,
    /// Committed to differentiate, and leaving the layer it divided in.
    Differentiating,
    /// Started its own programmed death, and breaking up.
    Apoptotic,
}

/// What happened in one lattice step.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StepReport {
    /// Index of the step just taken, counting from zero.
    pub mcs: u64,
    /// Relaxation sweeps the field took to settle at the start of this step.
    pub relax_steps: usize,
    /// Mean oxygen the first cell read at the start of this step.
    pub first_cell_oxygen: f64,
    /// Level of the lethal stress held at its source during this step: the
    /// exposure's intensity while it is on, zero otherwise or without one.
    pub exposure: f64,
    /// The tissue touched the edge of the box during this step.
    ///
    /// The bath is defined by what it can reach from the edge, so from here on
    /// the box decides what the tissue sees and everything after would be an
    /// artefact of its size. Division was skipped for this step, and a caller
    /// should stop.
    pub reached_edge: bool,
}

/// A lattice, the network inside each of its cells, and the oxygen around them,
/// stepped together.
pub struct Simulation<S: SimState> {
    scenario: Scenario,
    lattice: CpmLattice,
    /// One network per cell, or `None` once the cell is gone. Identifier `i + 1`
    /// names the cell whose slot is `i`, which holds as long as a daughter slot
    /// is pushed exactly when a daughter cell is created.
    networks: Vec<Option<Network>>,
    /// The clone of the cell in each slot, kept even once the cell is gone.
    lineage: Vec<CloneId>,
    /// Every clone ever founded, by identifier.
    clones: Vec<CloneRecord>,
    /// The network every first cell started with, unmutated, kept as the
    /// yardstick a mutant's response is measured against.
    pristine: Network,
    oxygen: ScalarField,
    /// The lethal stress, when the scenario exposes the tissue to one.
    stress: Option<ScalarField>,
    /// Level the stress field was last solved for, so that a change of level
    /// is solved from scratch, tightly, like the oxygen at the start.
    stress_level: f64,
    rng: StdRng,
    first: CellId,
    /// Index of the next step to take.
    next_mcs: u64,
    born: u64,
    gone: u64,
    _state: PhantomData<S>,
}

impl Simulation<Configured> {
    /// Seed the tissue and compile the hypothesis every cell will carry.
    ///
    /// Where the oxygen comes from decides where the first cell goes: at the
    /// centre of a bath, or against the wall of a vessel.
    ///
    /// # Errors
    ///
    /// [`CellweaveError::Config`] when the grid dimensions are invalid, the
    /// network does not compile, the mutation rate is not a probability, or a
    /// sheet is asked for with oxygen consumed or cells too small to divide.
    pub fn new(scenario: Scenario, spec: &NetworkSpec) -> Result<Self, CellweaveError> {
        if !(0.0..=1.0).contains(&scenario.mutation_rate) {
            return Err(CellweaveError::Config(format!(
                "the mutation rate is {}; it is a probability per division, between 0 and 1",
                scenario.mutation_rate
            )));
        }
        let mut lattice = CpmLattice::new(scenario.width, scenario.height, scenario.temperature)
            .map_err(|e| CellweaveError::Config(e.to_string()))?;
        let cx = scenario.width / 2;
        let cy = scenario.height / 2;

        if let OxygenSource::Vessel { radius } = scenario.oxygen_source {
            lattice.add_obstacle_disc(Pos2 { x: cx, y: cy }, radius);
        }
        let first = match scenario.seeding {
            Seeding::Single => match scenario.oxygen_source {
                OxygenSource::Bath => lattice.add_cell(CellKind::Tumor, Pos2 { x: cx, y: cy }, 8),
                OxygenSource::Vessel { radius } => lattice.add_cell(
                    CellKind::Tumor,
                    Pos2 {
                        x: cx + radius + 8,
                        y: cy,
                    },
                    8,
                ),
            },
            Seeding::Sheet { side } => {
                // A sheet leaves no medium for a bath to reach, and a field
                // that nothing supplies cannot be consumed without starving
                // everything; it runs with the field off.
                if scenario.oxygen_uptake != 0.0 {
                    return Err(CellweaveError::Config(
                        "a sheet fills the box and leaves no medium to supply oxygen; \
                         set oxygen_uptake = 0 to run it with the field off"
                            .into(),
                    ));
                }
                if side < 4 {
                    return Err(CellweaveError::Config(format!(
                        "a sheet of cells {side} pixels across is too fine for any of them to divide"
                    )));
                }
                // Every cell of a sheet is measured against the same reference
                // size, side by side. A cell clipped by the edge of the box
                // would read as crowded from the first step and be the first
                // to leave, so the box has to hold a whole number of them.
                if !scenario.width.is_multiple_of(side) || !scenario.height.is_multiple_of(side) {
                    return Err(CellweaveError::Config(format!(
                        "a sheet of cells {side} pixels across does not tile a {}x{} box; \
                         make both sides a multiple of {side}",
                        scenario.width, scenario.height
                    )));
                }
                // A vessel would leave blocks with no pixel at all, cells that
                // exist only to be counted as born and gone.
                if let OxygenSource::Vessel { .. } = scenario.oxygen_source {
                    return Err(CellweaveError::Config(
                        "a sheet is a layer of normal tissue with nothing inside it; \
                         it cannot be laid around a vessel"
                            .into(),
                    ));
                }
                let mut first = None;
                for y in (0..scenario.height).step_by(side as usize) {
                    for x in (0..scenario.width).step_by(side as usize) {
                        let id = lattice.add_block(CellKind::Normal, Pos2 { x, y }, side, side);
                        first.get_or_insert(id);
                    }
                }
                first.unwrap_or(CellId(1))
            }
        };
        let newborn_volume = match scenario.seeding {
            Seeding::Single => lattice.volume(first).map_or(0, |v| v.0),
            Seeding::Sheet { side } => side.checked_mul(side).ok_or_else(|| {
                CellweaveError::Config(format!(
                    "a sheet of cells {side} pixels across is too large"
                ))
            })?,
        };
        let seed_network = Network::from_spec(spec, newborn_volume)?;

        // A sheet labels every cell as the founder of its own clone, the way a
        // lineage tracing experiment marks every cell at the start and follows
        // what each one's descendants become.
        let cells = lattice.cell_count();
        let networks = vec![Some(seed_network.clone()); cells];
        let lineage = (0..cells as u32).map(CloneId).collect();
        let clones = (0..cells as u32)
            .map(|c| CloneRecord {
                id: CloneId(c),
                ..CloneRecord::founder()
            })
            .collect();

        // With nothing consuming, the field is not solved: everywhere sits at
        // the source's level, and starts there.
        let initial = match scenario.oxygen_source {
            _ if scenario.oxygen_uptake == 0.0 => 1.0,
            OxygenSource::Bath => 1.0,
            OxygenSource::Vessel { .. } => 0.0,
        };
        let oxygen = ScalarField::new(
            scenario.width,
            scenario.height,
            OXYGEN_DIFFUSION,
            0.0,
            initial,
        )?;
        let stress = match scenario.exposure {
            None => None,
            Some(e) => {
                // A field nothing takes up never settles, and fills the box
                // to the source's level whatever it is.
                if !(e.intensity.is_finite() && e.intensity >= 0.0) {
                    return Err(CellweaveError::Config(format!(
                        "the exposure's intensity is {}; it must be zero or more",
                        e.intensity
                    )));
                }
                if !(e.uptake.is_finite() && e.uptake > 0.0) {
                    return Err(CellweaveError::Config(format!(
                        "the exposure's uptake is {}; the cells must take some of it up, \
                         or it fills the box whatever its source",
                        e.uptake
                    )));
                }
                if let cellweave_core::Schedule::Periodic { on: 0, .. } = e.schedule {
                    return Err(CellweaveError::Config(
                        "a periodic exposure that is on for zero steps is never on".into(),
                    ));
                }
                if scenario.oxygen_uptake == 0.0 {
                    return Err(CellweaveError::Config(
                        "an exposure comes in where oxygen does, and this scenario has the \
                         field off; it has no source to come in from"
                            .into(),
                    ));
                }
                Some(ScalarField::new(
                    scenario.width,
                    scenario.height,
                    STRESS_DIFFUSION,
                    0.0,
                    0.0,
                )?)
            }
        };

        Ok(Self {
            rng: StdRng::seed_from_u64(scenario.seed),
            scenario,
            lattice,
            networks,
            lineage,
            clones,
            pristine: seed_network,
            oxygen,
            stress,
            stress_level: 0.0,
            first,
            next_mcs: 0,
            born: cells as u64,
            gone: 0,
            _state: PhantomData,
        })
    }

    /// Solve the field from scratch, tightly, so that every step after only
    /// ever has to nudge an equilibrium. Returns the running simulation and the
    /// number of sweeps solving took.
    ///
    /// # Errors
    ///
    /// [`CellweaveError::Simulation`] when the field does not settle.
    pub fn start(mut self) -> Result<(Simulation<Running>, usize), CellweaveError> {
        let sweeps = if self.field_is_off() {
            0
        } else {
            self.declare_sources();
            self.oxygen
                .relax(
                    TimeStep(1.0),
                    OXYGEN_WARMUP_TOLERANCE,
                    OXYGEN_WARMUP_MAX_STEPS,
                )
                .map_err(|e| {
                    CellweaveError::Simulation(format!("settling the initial field: {e}"))
                })?
        };
        let running = Simulation {
            scenario: self.scenario,
            lattice: self.lattice,
            networks: self.networks,
            lineage: self.lineage,
            clones: self.clones,
            pristine: self.pristine,
            oxygen: self.oxygen,
            stress: self.stress,
            stress_level: self.stress_level,
            rng: self.rng,
            first: self.first,
            next_mcs: self.next_mcs,
            born: self.born,
            gone: self.gone,
            _state: PhantomData,
        };
        Ok((running, sweeps))
    }
}

impl Simulation<Running> {
    /// Take one lattice step: settle the field, let every cell read it and run
    /// its network, decide fates, move pixels, and divide what is ready.
    ///
    /// # Errors
    ///
    /// [`CellweaveError::Simulation`] when the field does not settle.
    pub fn step(&mut self) -> Result<StepReport, CellweaveError> {
        let mcs = self.next_mcs;
        self.next_mcs += 1;

        // 1. Downward coupling: the field settles against the current layout of
        //    the cells, and each cell reads the oxygen where it actually sits.
        let relax_steps = if self.field_is_off() {
            0
        } else {
            self.declare_sources();
            self.oxygen
                .relax(
                    TimeStep(1.0),
                    OXYGEN_TOLERANCE * self.scenario.oxygen_uptake,
                    OXYGEN_MAX_STEPS,
                )
                .map_err(|e| CellweaveError::Simulation(format!("at MCS {mcs}: {e}")))?
        };
        let local_o2 = mean_per_cell(&self.oxygen, &self.lattice);
        let exposure = self.settle_stress(mcs)?;
        let local_stress = self
            .stress
            .as_ref()
            .map(|field| mean_per_cell(field, &self.lattice));

        // 2. Advance every living network by NETWORK_TIME_PER_MCS units of its
        //    time. A dead cell's network is frozen where it stopped, and so is
        //    that of a cell leaving the layer or breaking up: it has made its
        //    last decision.
        for (i, slot) in self.networks.iter_mut().enumerate() {
            let Some(net) = slot else { continue };
            let id = cell_of(i);
            if matches!(
                self.lattice.cell_kind(id),
                Some(CellKind::Necrotic | CellKind::Differentiated | CellKind::Apoptotic)
            ) {
                continue;
            }
            let inputs = CellInputs {
                o2: Concentration(local_o2[i + 1]),
                egf: Concentration(0.2),
                vegf: Concentration(0.0),
                crowding: Concentration(crowding(&self.lattice, id, net)),
                stress: Concentration(local_stress.as_ref().map_or(0.0, |s| s[i + 1])),
            };
            for _ in 0..NETWORK_TIME_PER_MCS {
                net.step(TimeStep(1.0), &inputs);
            }
        }

        // 3. Fate and upward coupling. A necrotic cell gives up pixels until
        //    none is left, and is then gone for good, and so does a cell leaving
        //    the layer or breaking up in apoptosis, faster; a cell whose
        //    survival just collapsed becomes necrotic; a cell whose network
        //    carries an apoptosis node may start apoptosis, and one whose network
        //    carries a differentiation node may commit to differentiate, each
        //    with that node's level as its chance this step; every other cell's
        //    network decides how large it tries to be.
        for (i, slot) in self.networks.iter_mut().enumerate() {
            let id = cell_of(i);
            let Some(net) = slot else { continue };
            let kind = self.lattice.cell_kind(id);
            match kind {
                Some(
                    leaving @ (CellKind::Necrotic | CellKind::Differentiated | CellKind::Apoptotic),
                ) => {
                    let volume = self.lattice.volume(id).map_or(0, |v| v.0);
                    if volume == 0 {
                        *slot = None;
                        self.gone += 1;
                    } else {
                        let rate = match leaving {
                            CellKind::Necrotic => RESORPTION_RATE,
                            CellKind::Apoptotic => FRAGMENTATION_RATE,
                            _ => SHEDDING_RATE,
                        };
                        let target = Volume(volume.saturating_sub(rate));
                        self.lattice
                            .set_volume_target(id, target, net.mechanics().lambda_volume);
                    }
                }
                Some(_) if net.survival() < DEATH_THRESHOLD => {
                    self.lattice.set_kind(id, CellKind::Necrotic);
                }
                // Each drawn only for a network that has the role, so a
                // network without it takes exactly the draws it took before
                // the role existed.
                Some(_)
                    if net
                        .apoptosis()
                        .is_some_and(|chance| self.rng.random::<f64>() < chance) =>
                {
                    self.lattice.set_kind(id, CellKind::Apoptotic);
                }
                Some(_)
                    if net
                        .differentiation()
                        .is_some_and(|chance| self.rng.random::<f64>() < chance) =>
                {
                    self.lattice.set_kind(id, CellKind::Differentiated);
                }
                Some(_) => {
                    let m = net.mechanics();
                    self.lattice
                        .set_volume_target(id, m.target_volume, m.lambda_volume);
                }
                None => {}
            }
        }

        // 4. Advance CPM.
        self.lattice.monte_carlo_step(&mut self.rng);

        let report = StepReport {
            mcs,
            relax_steps,
            first_cell_oxygen: local_o2.get(self.first.0 as usize).copied().unwrap_or(0.0),
            exposure,
            // A sheet touches the edge from the start; the edge only means
            // something for a tissue grown from one cell inside a supply.
            reached_edge: self.scenario.seeding == Seeding::Single && self.touches_the_edge(),
        };
        if report.reached_edge {
            return Ok(report);
        }

        // 5. Cells that have grown enough, and are not too short of oxygen to
        //    cycle, split in two.
        let every = match self.scenario.seeding {
            Seeding::Single => DIVISION_CHECK_EVERY,
            Seeding::Sheet { .. } => 1,
        };
        if mcs.is_multiple_of(every) && self.living() < MAX_CELLS {
            self.born += self.divide_ready_cells(mcs);
        }

        Ok(report)
    }

    /// Split every cell, tumour or normal, that has grown past the division
    /// threshold and is still cycling, and say how many daughters were born.
    ///
    /// Hypoxia arrests the cell cycle long before it kills: a cell whose hypoxia
    /// response is high is quiescent, alive and holding its place, but not
    /// dividing. Without that middle state nothing bounds growth short of death.
    ///
    /// A daughter inherits a copy of its parent's network, so its protein state
    /// and its parameters, which is what makes the population a lineage rather
    /// than a crowd of identical newborns. With the scenario's mutation rate, the
    /// copy has one rate changed, and the daughter founds a clone.
    fn divide_ready_cells(&mut self, mcs: u64) -> u64 {
        let ready: Vec<usize> = self
            .networks
            .iter()
            .enumerate()
            .filter_map(|(i, slot)| slot.as_ref().map(|net| (i, net)))
            .filter(|(i, net)| {
                let id = cell_of(*i);
                let threshold = DIVISION_RATIO * f64::from(net.reference_volume());
                let big_enough = self
                    .lattice
                    .volume(id)
                    .is_some_and(|v| f64::from(v.0) >= threshold);
                big_enough && self.fate(id) == Some(Fate::Cycling)
            })
            .map(|(i, _)| i)
            .collect();

        let mut born = 0;
        let living = self.living();
        for i in ready {
            if living + born as usize >= MAX_CELLS {
                break;
            }
            let angle = self.rng.random::<f64>() * std::f64::consts::PI;
            if self.lattice.divide(cell_of(i), angle).is_some() {
                let mut daughter = self.networks[i].clone();
                let mut clone = self.lineage[i];
                // The rate is tested before anything is drawn, so that a run
                // without mutation takes exactly the draws it always took.
                let rate = self.scenario.mutation_rate;
                if rate > 0.0 && self.rng.random::<f64>() < rate {
                    let (which, how) = (self.rng.random::<f64>(), self.rng.random::<f64>());
                    let founded = CloneId(self.clones.len() as u32);
                    let suppressed = self.scenario.suppress.contains(&founded);
                    // A suppressed mutation is made on a copy, to know what it
                    // would have been, and the copy is thrown away.
                    let made = if suppressed {
                        daughter.clone().and_then(|mut n| n.mutate(which, how))
                    } else {
                        daughter.as_mut().and_then(|n| n.mutate(which, how))
                    };
                    if let Some(mutation) = made {
                        self.clones.push(CloneRecord {
                            id: founded,
                            parent: Some(clone),
                            founded_at: mcs,
                            mutation: Some(mutation),
                            suppressed,
                        });
                        clone = founded;
                    }
                }
                self.networks.push(daughter);
                self.lineage.push(clone);
                born += 1;
            }
        }

        // Kept as a hard assertion rather than a comment: were the two to drift,
        // a network would quietly drive the wrong cell and every result after
        // that would be wrong without a single test failing.
        assert_eq!(
            self.networks.len(),
            self.lattice.cell_count(),
            "one slot per cell: the lattice and the network list have diverged"
        );
        assert_eq!(
            self.networks.len(),
            self.lineage.len(),
            "one clone per cell: the network list and the lineage have diverged"
        );
        born
    }

    /// Every cell still holding pixels, with the oxygen it sits in, labelled
    /// with the step last taken.
    ///
    /// Daughters born this step have not read the field yet; dividing only
    /// relabels pixels, so reading the field again against the current layout
    /// is exact. A cell that is gone has nothing left to record.
    pub fn snapshot(&self) -> SimSnapshot {
        let local_o2 = mean_per_cell(&self.oxygen, &self.lattice);
        let cells = self
            .lattice
            .cell_records()
            .zip(self.networks.iter())
            .filter_map(|((id, rec), slot)| {
                let net = slot.as_ref()?;
                let clone = self.clone_of(id)?;
                (rec.volume.0 > 0).then(|| CellSnapshot {
                    id,
                    kind: rec.kind,
                    volume: rec.volume,
                    proteins: net.state().clone(),
                    oxygen: Concentration(local_o2[id.0 as usize]),
                    clone,
                })
            })
            .collect();
        SimSnapshot {
            mcs: self.last_mcs(),
            cells,
        }
    }
}

impl<S: SimState> Simulation<S> {
    /// Solve the stress field for step `mcs`, and say what level its source
    /// held. Without an exposure there is no field and the level is zero.
    ///
    /// While exposure is off the field is zero everywhere: nothing supplies
    /// it, and everything that takes it up has had it. When the level changes,
    /// switching on or off, the field is solved again from scratch and
    /// tightly, for the reason the oxygen is at the start: a per-step
    /// tolerance lets a slow front creep on unnoticed.
    fn settle_stress(&mut self, mcs: u64) -> Result<f64, CellweaveError> {
        let (Some(exposure), Some(field)) = (self.scenario.exposure, self.stress.as_mut()) else {
            return Ok(0.0);
        };
        let level = exposure.level_at(mcs);
        if level == 0.0 {
            field.fill(0.0);
            self.stress_level = 0.0;
            return Ok(0.0);
        }
        apply_uptake(field, &self.lattice, exposure.uptake);
        field.release_all();
        match self.scenario.oxygen_source {
            OxygenSource::Bath => hold_medium_bath(field, &self.lattice, level),
            OxygenSource::Vessel { .. } => {
                apply_background_uptake(field, &self.lattice, exposure.uptake);
                hold_obstacles_at(field, &self.lattice, level);
            }
        }
        let (tolerance, max_sweeps) = if level == self.stress_level {
            (OXYGEN_TOLERANCE * exposure.uptake, OXYGEN_MAX_STEPS)
        } else {
            (OXYGEN_WARMUP_TOLERANCE, OXYGEN_WARMUP_MAX_STEPS)
        };
        field
            .relax(TimeStep(1.0), tolerance, max_sweeps)
            .map_err(|e| CellweaveError::Simulation(format!("stress at MCS {mcs}: {e}")))?;
        self.stress_level = level;
        Ok(level)
    }

    /// Whether the field is left unsolved: nothing consumes, so everywhere
    /// sits at the level it started at.
    fn field_is_off(&self) -> bool {
        self.scenario.oxygen_uptake == 0.0
    }

    /// Tell the field who consumes and who supplies, from the current layout.
    ///
    /// Around a bath, living cells consume and the medium the bath can reach is
    /// held full. Around a vessel there is no empty liquid but stroma, so the
    /// medium consumes too, and the vessel wall is what is held.
    fn declare_sources(&mut self) {
        let uptake = self.scenario.oxygen_uptake;
        apply_uptake(&mut self.oxygen, &self.lattice, uptake);
        self.oxygen.release_all();
        match self.scenario.oxygen_source {
            OxygenSource::Bath => hold_medium_bath(&mut self.oxygen, &self.lattice, 1.0),
            OxygenSource::Vessel { .. } => {
                apply_background_uptake(&mut self.oxygen, &self.lattice, uptake);
                hold_obstacles_at(&mut self.oxygen, &self.lattice, 1.0);
            }
        }
    }

    /// Whether any cell occupies a pixel on the edge of the box.
    fn touches_the_edge(&self) -> bool {
        let lattice = &self.lattice;
        let (w, h) = lattice.dims();
        let is_cell = |x: u32, y: u32| {
            let id = lattice.occupant(x, y);
            id != 0 && id != OBSTACLE
        };
        (0..w).any(|x| is_cell(x, 0) || is_cell(x, h - 1))
            || (0..h).any(|y| is_cell(0, y) || is_cell(w - 1, y))
    }

    /// Index of the step last taken, which labels what is read off the
    /// simulation now. Zero before the first step as well as after it.
    pub(crate) fn last_mcs(&self) -> u64 {
        self.next_mcs.saturating_sub(1)
    }

    /// What this simulation was asked to be.
    pub fn scenario(&self) -> &Scenario {
        &self.scenario
    }

    /// The lattice: which cell holds each pixel, and each cell's record.
    pub fn lattice(&self) -> &CpmLattice {
        &self.lattice
    }

    /// The lethal stress field, as it was last settled, or `None` when the
    /// scenario exposes the tissue to nothing.
    pub fn stress(&self) -> Option<&ScalarField> {
        self.stress.as_ref()
    }

    /// The oxygen field, as it was last settled.
    pub fn oxygen(&self) -> &ScalarField {
        &self.oxygen
    }

    /// The network inside a cell, or `None` for the medium, an obstacle, or a
    /// cell that is gone.
    pub fn network(&self, id: CellId) -> Option<&Network> {
        let slot = (id.0 as usize).checked_sub(1)?;
        self.networks.get(slot)?.as_ref()
    }

    /// The clone a cell belongs to, or `None` for the medium or an obstacle. A
    /// cell that is gone keeps the clone it had.
    pub fn clone_of(&self, id: CellId) -> Option<CloneId> {
        let slot = (id.0 as usize).checked_sub(1)?;
        self.lineage.get(slot).copied()
    }

    /// The network the first cells started with, before any mutation.
    pub fn pristine(&self) -> &Network {
        &self.pristine
    }

    /// Every clone founded so far, the founding clone first, indexed by
    /// identifier.
    pub fn clones(&self) -> &[CloneRecord] {
        &self.clones
    }

    /// What a cell is doing, or `None` for the medium, an obstacle, or a cell
    /// that is gone.
    ///
    /// A tumour or normal cell cycles while its hypoxia response is below the
    /// quiescence threshold; a cell committed to differentiate is leaving, an
    /// apoptotic one is breaking up, and a necrotic one is dead.
    pub fn fate(&self, id: CellId) -> Option<Fate> {
        let net = self.network(id)?;
        Some(match self.lattice.cell_kind(id)? {
            CellKind::Necrotic => Fate::Necrotic,
            CellKind::Differentiated => Fate::Differentiating,
            CellKind::Apoptotic => Fate::Apoptotic,
            CellKind::Tumor | CellKind::Normal if net.hypoxia_response() < QUIESCENCE_THRESHOLD => {
                Fate::Cycling
            }
            _ => Fate::Quiescent,
        })
    }

    /// Cells alive and in the layer: neither necrotic, nor leaving, nor
    /// breaking up, and still holding pixels.
    pub fn living(&self) -> usize {
        self.lattice
            .cell_records()
            .filter(|(_, r)| {
                !matches!(
                    r.kind,
                    CellKind::Necrotic | CellKind::Differentiated | CellKind::Apoptotic
                ) && r.volume.0 > 0
            })
            .count()
    }

    /// Necrotic cells not yet resorbed.
    pub fn necrotic(&self) -> usize {
        self.lattice
            .cell_records()
            .filter(|(_, r)| r.kind == CellKind::Necrotic && r.volume.0 > 0)
            .count()
    }

    /// Cells ever born, the first one included.
    pub fn born(&self) -> u64 {
        self.born
    }

    /// Cells that have left the tissue for good: resorbed to nothing after
    /// dying, or shed after differentiating.
    pub fn gone(&self) -> u64 {
        self.gone
    }
}

/// How squeezed a cell is: how far below its lineage's reference volume it
/// sits, as a fraction of it. Zero for a cell at or above that size.
///
/// What a cell can sense of its neighbours' pressure is its own spread: a cell
/// pressed small by the tissue around it has less room than it was born to,
/// and that is the signal contact inhibition in an epithelium runs on.
fn crowding(lattice: &CpmLattice, id: CellId, net: &Network) -> f64 {
    let volume = lattice.volume(id).map_or(0, |v| v.0);
    let reference = f64::from(net.reference_volume().max(1));
    (1.0 - f64::from(volume) / reference).clamp(0.0, 1.0)
}

/// The cell whose slot is `i`. Identifier 0 is the medium.
fn cell_of(i: usize) -> CellId {
    CellId(i as u32 + 1)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn spec() -> NetworkSpec {
        toml::from_str(include_str!("../../../examples/networks/ras-erk.toml"))
            .expect("the shipped network parses")
    }

    fn scenario(oxygen_source: OxygenSource) -> Scenario {
        Scenario {
            width: 80,
            height: 80,
            temperature: 10.0,
            oxygen_source,
            oxygen_uptake: 0.0003,
            seed: DEFAULT_SEED,
            mutation_rate: 0.0,
            suppress: Vec::new(),
            seeding: Seeding::Single,
            exposure: None,
        }
    }

    #[test]
    fn one_seed_gives_one_run() {
        let run = || {
            let (mut sim, _) = Simulation::new(scenario(OxygenSource::Bath), &spec())
                .unwrap()
                .start()
                .unwrap();
            for _ in 0..60 {
                sim.step().unwrap();
            }
            // Debug prints every float to the shortest string that reads back
            // exactly, so equal strings are equal runs.
            format!("{:?}", sim.snapshot())
        };
        assert_eq!(run(), run());
    }

    #[test]
    fn the_first_cell_sits_against_the_vessel() {
        let sim = Simulation::new(scenario(OxygenSource::Vessel { radius: 6 }), &spec()).unwrap();
        // Centre 40, radius 6, seeded 8 further out.
        let id = sim.lattice().occupant(54, 40);
        assert_eq!(id, 1, "the first cell is where the vessel puts it");
        assert_eq!(sim.lattice().occupant(40, 40), OBSTACLE);
    }

    #[test]
    fn a_network_is_found_by_its_cell_and_not_by_the_medium() {
        let sim = Simulation::new(scenario(OxygenSource::Bath), &spec()).unwrap();
        assert!(sim.network(CellId(1)).is_some());
        assert!(
            sim.network(CellId(0)).is_none(),
            "the medium has no network"
        );
        assert!(sim.network(CellId(2)).is_none(), "no second cell yet");
    }

    #[test]
    fn a_fresh_cell_in_full_oxygen_is_cycling() {
        let (sim, _) = Simulation::new(scenario(OxygenSource::Bath), &spec())
            .unwrap()
            .start()
            .unwrap();
        assert_eq!(sim.fate(CellId(1)), Some(Fate::Cycling));
        assert_eq!(sim.fate(CellId(0)), None, "the medium has no fate");
    }

    #[test]
    fn steps_are_counted_from_zero_and_label_the_snapshot() {
        let (mut sim, _) = Simulation::new(scenario(OxygenSource::Bath), &spec())
            .unwrap()
            .start()
            .unwrap();
        assert_eq!(sim.step().unwrap().mcs, 0);
        assert_eq!(sim.step().unwrap().mcs, 1);
        assert_eq!(sim.snapshot().mcs, 1);
    }

    fn evolving(rate: f64) -> Simulation<Running> {
        let scenario = Scenario {
            mutation_rate: rate,
            ..scenario(OxygenSource::Bath)
        };
        let (mut sim, _) = Simulation::new(scenario, &spec()).unwrap().start().unwrap();
        for _ in 0..80 {
            sim.step().unwrap();
        }
        sim
    }

    #[test]
    fn without_mutation_every_cell_is_the_founder() {
        let sim = evolving(0.0);
        assert!(sim.born() > 1, "the population grew");
        assert_eq!(sim.clones().len(), 1);
        assert!(sim.snapshot().cells.iter().all(|c| c.clone == CloneId(0)));
    }

    #[test]
    fn when_every_daughter_mutates_every_daughter_founds_a_clone() {
        let sim = evolving(1.0);
        assert_eq!(sim.clones().len() as u64, sim.born());
        for (i, clone) in sim.clones().iter().enumerate().skip(1) {
            assert_eq!(clone.id, CloneId(i as u32));
            let parent = clone.parent.unwrap();
            assert!(parent < clone.id, "a clone arises in an older one");
            assert!(clone.mutation.is_some());
        }
        assert_eq!(
            format!("{:?}", sim.snapshot()),
            format!("{:?}", evolving(1.0).snapshot()),
            "one seed, one evolution"
        );
    }

    #[test]
    fn a_mutation_rate_is_a_probability() {
        for rate in [-0.1, 1.5, f64::NAN] {
            let scenario = Scenario {
                mutation_rate: rate,
                ..scenario(OxygenSource::Bath)
            };
            assert!(Simulation::new(scenario, &spec()).is_err(), "{rate}");
        }
    }

    /// Every step's measurements, and who holds each pixel at the end.
    fn history(
        scenario: Scenario,
        steps: u32,
    ) -> (Vec<crate::Measurements>, Vec<u32>, Vec<CloneRecord>) {
        let (mut sim, _) = Simulation::new(scenario, &spec()).unwrap().start().unwrap();
        let mut measured = Vec::new();
        for _ in 0..steps {
            sim.step().unwrap();
            measured.push(sim.measure());
        }
        let (w, h) = sim.lattice().dims();
        let grid = (0..h)
            .flat_map(|y| (0..w).map(move |x| (x, y)))
            .map(|(x, y)| sim.lattice().occupant(x, y))
            .collect();
        (measured, grid, sim.clones().to_vec())
    }

    fn mutating() -> Scenario {
        Scenario {
            mutation_rate: 1.0,
            ..scenario(OxygenSource::Bath)
        }
    }

    #[test]
    fn a_suppressed_mutation_still_founds_its_clone() {
        let (_, _, real) = history(mutating(), 100);
        let target = real[3].id;
        let (_, _, without) = history(
            Scenario {
                suppress: vec![target],
                ..mutating()
            },
            100,
        );
        let kept = &without[target.0 as usize];
        assert!(kept.suppressed);
        assert_eq!(kept.parent, real[3].parent);
        assert_eq!(kept.founded_at, real[3].founded_at);
        assert_eq!(kept.mutation, real[3].mutation, "what it would have been");
        assert!(without.iter().filter(|c| c.suppressed).count() == 1);
    }

    #[test]
    fn two_histories_are_one_until_the_suppressed_mutation() {
        let (real, _, clones) = history(mutating(), 100);
        let target = clones.last().unwrap();
        let (without, _, _) = history(
            Scenario {
                suppress: vec![target.id],
                ..mutating()
            },
            100,
        );
        let before = target.founded_at as usize;
        assert_eq!(real[..=before], without[..=before]);
    }

    #[test]
    fn a_neutral_mutation_leaves_the_history_identical_to_the_last_bit() {
        // VEGF is read by nothing in the shipped network, so a change to how
        // it is made or lost cannot reach a fate.
        let (real, grid, clones) = history(mutating(), 100);
        let neutral = clones
            .iter()
            .find(|c| c.mutation.as_ref().is_some_and(|m| m.node == "vegf"))
            .expect("one of the clones changed VEGF");
        let (without, other_grid, _) = history(
            Scenario {
                suppress: vec![neutral.id],
                ..mutating()
            },
            100,
        );
        assert_eq!(real, without);
        assert_eq!(grid, other_grid);
    }
}
