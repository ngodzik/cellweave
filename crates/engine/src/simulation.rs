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

use std::collections::{BTreeMap, BTreeSet};
use std::marker::PhantomData;

use cellweave_core::{
    Angiogenesis, CellId, CellInputs, CellKind, CellSnapshot, CellweaveError, CloneId,
    Concentration, Cycle, DiffusionField, Exposure, Immunity, NecroticFactor, NetworkSpec,
    OxygenSource, Perfusion, Pos2, Pos3, Seeding, SignalingNetwork, SimSnapshot, StillMedium,
    TimeStep, Tracing, Volume, Waste,
};
use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};

use crate::coupling::{
    apply_background_uptake, apply_production, apply_uptake, hold_medium_bath, hold_vessels_at,
    hold_vessels_by, mean_per_cell,
};
use crate::cpm::{CpmLattice, OBSTACLE};
use crate::evolution::CloneRecord;
use crate::network::Network;
use crate::pde::ScalarField;
use crate::sisters::{Lives, SisterFates};
use crate::vessels::{Capillary, Tip, Vessels};

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

/// Refuse to go past this many living cells unless a scenario says
/// otherwise, as a guard against runaway growth. Counted among the living
/// rather than among every cell ever born, since a tissue that turns over
/// keeps its size while its births add up.
pub const DEFAULT_MAX_CELLS: usize = 2000;

/// Radius of the first cell in the plane, in pixels: a disc of about two
/// hundred pixels, which every later cell is measured against.
const FIRST_CELL_RADIUS: u32 = 8;

/// Radius of the first cell in a volume, in voxels: a ball of 123 voxels.
///
/// A volume cannot be cut as finely as a plane. A cell about 16 pixels
/// across in the plane would be some two thousand voxels, and a spheroid as
/// large as measured ones tens of thousands of cells, a lattice of billions.
/// At three, a cell is about seven voxels across, a voxel some 2.5 to 3 µm,
/// and the measured spheroids fit in a box of a few hundred on a side.
const FIRST_CELL_RADIUS_3D: u32 = 3;

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

/// How many times a step the edge of dead tissue is drawn from while its
/// debris owes pixels, each time laying bare the dead pixels behind those
/// taken.
///
/// In Greenspan's balance the edge of a core moves in at $\lambda \eta R / 3$
/// a step: at a loss of 0.16 a step, about three voxels a step for a core of
/// radius 60. Eight leave room, and a core that owes more carries it.
const EDGE_PASSES: usize = 8;

/// The lethal stress diffuses like oxygen. How far it reaches into the tissue
/// is set by how much the cells take up, which the exposure states, so one
/// diffusion coefficient for both keeps the comparison between the two fields
/// a comparison of their uptakes alone.
const STRESS_DIFFUSION: f64 = OXYGEN_DIFFUSION;

/// The waste diffuses like oxygen. In Greenspan's model only its production
/// against its diffusion, and the level that arrests a cell, decide where
/// cells arrest, so one coefficient for both leaves the production the one
/// rate to set.
const WASTE_DIFFUSION: f64 = OXYGEN_DIFFUSION;

/// The necrotic factor diffuses as the waste does. Jiang and colleagues put
/// such a factor's diffusion tens of thousands of times below oxygen's, that
/// of a protein of 80 to 90 kilodaltons; but a field held at equilibrium, as
/// every field here is, depends on its production over its diffusion alone,
/// and the production is the rate a scenario sets.
const NECROTIC_FACTOR_DIFFUSION: f64 = WASTE_DIFFUSION;

/// VEGF diffuses like oxygen; how far it reaches is set by how fast it
/// decays, which the scenario's reach states: $\lambda = D / L^2$.
const VEGF_DIFFUSION: f64 = OXYGEN_DIFFUSION;

/// The VEGF field is relaxed until no pixel moves by more than this fraction
/// of what a cell secreting at a level of one adds per step, for the reason
/// the oxygen's tolerance is relative to its uptake.
const VEGF_TOLERANCE: f64 = 0.01;

/// Radius of an immune cell where it arrives, in pixels: a disc of 13
/// pixels, about a third of a tumour cell's area, a lymphocyte being small.
const IMMUNE_RADIUS: u32 = 2;

/// The size an immune cell keeps to, in pixels.
const IMMUNE_VOLUME: u32 = 13;

/// How firmly an immune cell keeps to its size: loosely, since a crawling
/// lymphocyte changes shape as it goes. A tumour cell's 50 would hold so
/// small a cell still, every pixel it gained or lost costing more than the
/// chemokine's pull could pay.
const IMMUNE_LAMBDA: f64 = 2.0;

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
    /// Grid depth in voxels: one for a tissue in the plane, more for one in
    /// a volume.
    pub depth: u32,
    /// The most living cells the run lets divide.
    pub max_cells: usize,
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
    /// Probability that a mutation changes the network's wiring, gaining or
    /// losing a term, rather than scaling a rate. Zero, the default, draws
    /// nothing for it, so a run without rewiring is the run it was before.
    pub rewiring: f64,
    /// Clones whose founding mutation is, in this run, drawn but not made.
    ///
    /// Everything else happens as it would have: the same numbers are drawn,
    /// the clone is founded with the same identifier and the same parent, and
    /// its cells are the cells they would have been. Only their network is
    /// their parent's, unchanged. Until that difference makes a difference,
    /// the run is the run without the suppression, bit for bit; a mutation
    /// that never makes one leaves it so to the end.
    pub suppress: Vec<CloneId>,
    /// Every mutation founded before this step is drawn but not made, as
    /// [`Scenario::suppress`] does for one clone: the tissue has the same
    /// luck and none of the mutant clones. Left out, nothing is suppressed
    /// this way.
    pub suppress_before: Option<u64>,
    /// Probability that a daughter is born with one rate of its network
    /// changed, founding a clone of its own. Zero keeps every daughter an
    /// exact copy, and draws nothing from the generator for it, so a run
    /// without mutation is the same run it was before mutation existed.
    pub mutation_rate: f64,
    /// An early tumour that arises in a sheet at one step, if any. Left out,
    /// the run is the one it always was.
    pub lesion: Option<Lesion>,
    /// Clones whose founding mutation is undone at one step, if any. Left
    /// out, nothing is undone.
    pub revert: Option<Reversion>,
    /// Vessels that grow towards the tissue that calls for them. Left out,
    /// the vessel stays as it was laid and nothing reads what cells secrete.
    pub angiogenesis: Option<Angiogenesis>,
    /// Immune cells that come for the tumour. Left out, none do.
    pub immunity: Option<Immunity>,
    /// Every cell marked at one step and the clones the marks make
    /// followed. Left out, nothing is marked.
    pub tracing: Option<Tracing>,
    /// The share of a dead cell's volume its debris loses per step, broken
    /// down into fluid that leaves the tissue, as the necrotic core of a
    /// spheroid loses mass (Greenspan 1972). The living tissue closes in
    /// over what goes, at the edge of the dead tissue, so that what is lost
    /// in the core is given up at the surface. Which pixels go, and to which
    /// neighbour, is drawn from the run's generator. Left out, a dead cell
    /// loses a pixel only when a neighbour takes it, which in a core of dead
    /// cells none does, and the run is the one it always was.
    pub debris_loss: Option<f64>,
    /// Monte Carlo sweeps of the lattice each step: how many times every
    /// pixel that can change is tried before the networks, the fields and
    /// the fates move on. One in every configuration calibrated so far; more
    /// lets a packed tissue's mechanics keep pace with what its growth and
    /// its debris loss ask of it (chapter 27).
    pub sweeps_per_step: u32,
    /// A cell cycle that takes time, each cell drawing its own duration at
    /// birth from the run's generator and growing through it; a tissue grown
    /// from one cell is then checked for division every step, each cell
    /// against its own clock. Left out, a cell asks at birth for all it will
    /// be and divides once it is big enough at a check, and the run is the
    /// one it always was.
    pub cycle: Option<Cycle>,
    /// How far two pixels count as touching, for the contact energy: 1
    /// across a face, as in every configuration before chapter 27; in a
    /// volume 3 across a face, an edge or a corner, so that what a boundary
    /// costs depends less on its direction ([`CpmLattice::set_contact_order`]).
    /// Copies are proposed across faces whatever it is.
    pub contact_order: u8,
    /// A waste the living release and the bath washes away, which each
    /// cell's network reads as its `waste` input. `None` leaves the field
    /// unsolved and every cell reading zero, and the run is the one it was.
    pub waste: Option<Waste>,
    /// A factor the dead release and the bath washes away, which each cell's
    /// network reads as its `necrotic_factor` input. `None` leaves the field
    /// unsolved and every cell reading zero, and the run is the one it was.
    pub necrotic_factor: Option<NecroticFactor>,
    /// A bath whose medium is still, so that what reaches the tissue's
    /// surface is what diffuses through it: the more the tissue consumes, the
    /// less. Only around a spheroid in a volume. `None` holds the medium at
    /// the source's level up to the tissue, as if stirred, and the run is the
    /// one it was.
    pub still_medium: Option<StillMedium>,
}

/// One immune cell in the tissue.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct ImmuneCell {
    /// Its slot, one less than its identifier.
    slot: usize,
    /// The step it arrived at.
    arrived: u64,
    /// How many tumour cells it can still kill.
    kills_left: u32,
    /// Whether it has died, exhausted or old, and is breaking up.
    dying: bool,
}

/// Clones whose founding mutation is undone, at one step, in every living
/// cell that carries it.
///
/// [`Scenario::suppress`] asks what would have happened had a mutation never
/// been made; this asks what happens if it stops acting, here and now. The
/// run is the run without the reversion, to the last bit, until the step it
/// comes at, and from then on differs only by those cells being as their
/// lineage would have been without that mutation. It is what an experiment
/// does when it switches a gene back on in the cells of one clone.
#[derive(Debug, Clone, PartialEq)]
pub struct Reversion {
    /// The clones whose founding mutation is undone, in them and in every
    /// clone that arose in them.
    pub clones: Vec<CloneId>,
    /// The step it is undone at, before anything else in that step.
    pub at: u64,
}

/// A patch of tissue whose cells, at one step, start carrying another
/// network: in a sheet, an early tumour arising inside normal tissue; in a
/// tumour, a subclone placed by hand.
///
/// What sets it apart is only what its network does, its cells growing,
/// dividing, leaving or dying as every other cell does. That is all a
/// lesion needs to be for the questions it is here to ask: whether its
/// neighbours outcompete it, or how much of what it does goes to them.
#[derive(Debug, Clone, PartialEq)]
pub struct Lesion {
    /// The network its cells carry from the step it appears.
    pub spec: NetworkSpec,
    /// The step it appears at, before anything else happens in that step.
    pub at: u64,
    /// Every living cell whose centre lies within this many pixels of the
    /// lesion's centre joins it.
    pub radius: f64,
    /// Where the lesion is centred, in pixels across and down the box; the
    /// centre of the box when left out.
    pub centre: Option<(f64, f64)>,
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
    /// What each cell a division made did next, and its sister.
    lives: Lives,
    /// The network a lesion's cells start carrying, compiled once, when the
    /// scenario has a lesion.
    lesion_network: Option<Network>,
    oxygen: ScalarField,
    /// The lethal stress, when the scenario exposes the tissue to one.
    stress: Option<ScalarField>,
    /// Level the stress field was last solved for, so that a change of level
    /// is solved from scratch, tightly, like the oxygen at the start.
    stress_level: f64,
    /// What the living secrete as a call for vessels, when vessels can grow.
    vegf: Option<ScalarField>,
    /// The waste the living release, when the scenario has one.
    waste: Option<ScalarField>,
    /// The factor the dead release, when the scenario has one.
    necrotic_factor: Option<ScalarField>,
    /// The sprouts growing, and where every one started.
    vessels: Vessels,
    rng: StdRng,
    first: CellId,
    /// Index of the next step to take.
    next_mcs: u64,
    born: u64,
    gone: u64,
    /// Cells a growing vessel took every pixel of.
    crushed: u64,
    /// The chemokine living tumour cells secrete, when immune cells come.
    chemokine: Option<ScalarField>,
    /// The immune cells in the tissue, dying ones included.
    immune: Vec<ImmuneCell>,
    /// Tumour cells immune cells have killed.
    killed: u64,
    /// The mark each slot's cell carries, if tracing has marked it or an
    /// ancestor: the slot of the cell marked at the tracing step.
    marks: Vec<Option<u32>>,
    /// With debris loss, the slots of dead cells holding to what they had
    /// when they died.
    kept_dead: BTreeSet<usize>,
    /// The step each dead cell's slot died at, while the cell is there; for
    /// a cell made necrotic outside the fate step, the first step it was
    /// seen dead.
    died_at: BTreeMap<usize, u64>,
    /// With debris loss, the pixels the dead tissue owes the fluid and has
    /// not yet given up.
    debris_owed: f64,
    /// The level the bath holds the medium at, next to the tissue: the
    /// source's, 1, unless the medium is still.
    bath_level: f64,
    /// With a cycle, the step from which each slot's cell may divide; zero
    /// for the cells the tissue starts with, whose cycle is behind them.
    due: Vec<u64>,
    /// With a cycle, the step each slot's cell was born at, its cycle
    /// began; zero for the cells the tissue starts with.
    began: Vec<u64>,
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
    /// network does not compile, the mutation rate is not a probability, a
    /// sheet is asked for with oxygen consumed or cells too small to divide,
    /// or tracing writes clones down every zero steps, or a cycle lasts less
    /// than a step or spreads by less than nothing, or the contact order is
    /// not 1, 2 or 3, or a waste or a necrotic factor is released at no rate
    /// or with no bath to wash it away.
    pub fn new(scenario: Scenario, spec: &NetworkSpec) -> Result<Self, CellweaveError> {
        if scenario.sweeps_per_step == 0 {
            return Err(CellweaveError::Config(
                "sweeps_per_step is 0; a step sweeps the lattice at least once".into(),
            ));
        }
        if !(1..=3).contains(&scenario.contact_order) {
            return Err(CellweaveError::Config(format!(
                "the contact order is {}; it is 1, 2 or 3",
                scenario.contact_order
            )));
        }
        if let Some(c) = scenario.cycle
            && !(c.median.is_finite() && c.median >= 1.0 && c.spread.is_finite() && c.spread >= 0.0)
        {
            return Err(CellweaveError::Config(format!(
                "the cycle's median is {} steps and its spread {}; a cycle lasts a step at \
                 least, and its spread is not below zero",
                c.median, c.spread
            )));
        }
        if !(0.0..=1.0).contains(&scenario.mutation_rate) {
            return Err(CellweaveError::Config(format!(
                "the mutation rate is {}; it is a probability per division, between 0 and 1",
                scenario.mutation_rate
            )));
        }
        if !(0.0..=1.0).contains(&scenario.rewiring) {
            return Err(CellweaveError::Config(format!(
                "rewiring is {}; it is the probability that a mutation rewires, between 0 and 1",
                scenario.rewiring
            )));
        }
        if let Some(loss) = scenario.debris_loss
            && !(loss.is_finite() && loss > 0.0 && loss <= 1.0)
        {
            return Err(CellweaveError::Config(format!(
                "debris loss is {loss}; it is the share of a dead cell lost per step, above \
                 zero and at most one"
            )));
        }
        if let Some(w) = scenario.waste {
            if !(w.production.is_finite() && w.production > 0.0) {
                return Err(CellweaveError::Config(format!(
                    "the waste's production is {}; the living release some, above zero",
                    w.production
                )));
            }
            if scenario.oxygen_source != OxygenSource::Bath {
                return Err(CellweaveError::Config(
                    "a waste is washed away by the bath, and this scenario has none".into(),
                ));
            }
        }
        if let Some(f) = scenario.necrotic_factor {
            if !(f.production.is_finite() && f.production > 0.0) {
                return Err(CellweaveError::Config(format!(
                    "the necrotic factor's production is {}; the dead release some, above zero",
                    f.production
                )));
            }
            if f.release_for == Some(0) {
                return Err(CellweaveError::Config(
                    "a necrotic factor released for 0 steps is never released; leave release_for \
                     out to release it for as long as a dead cell is there"
                        .into(),
                ));
            }
            if scenario.oxygen_source != OxygenSource::Bath {
                return Err(CellweaveError::Config(
                    "a necrotic factor is washed away by the bath, and this scenario has none"
                        .into(),
                ));
            }
        }
        if let Some(m) = scenario.still_medium {
            if !(m.diffusion.is_finite() && m.diffusion > 0.0) {
                return Err(CellweaveError::Config(format!(
                    "the still medium's diffusion is {}; oxygen diffuses through it, above zero",
                    m.diffusion
                )));
            }
            if scenario.oxygen_source != OxygenSource::Bath {
                return Err(CellweaveError::Config(
                    "a still medium is the bath's, and this scenario has none".into(),
                ));
            }
            if scenario.depth == 1 {
                return Err(CellweaveError::Config(
                    "a still medium is solved around a ball, and a plane has none; give the \
                     scenario a depth"
                        .into(),
                ));
            }
            // What the tissue releases would gather in a still medium as
            // oxygen is depleted in it, and only oxygen's level is solved.
            if scenario.waste.is_some() || scenario.necrotic_factor.is_some() {
                return Err(CellweaveError::Config(
                    "a still medium holds back what the tissue releases as well, and only \
                     oxygen is solved in it; leave out the waste and the necrotic factor"
                        .into(),
                ));
            }
        }
        if scenario.tracing.is_some_and(|t| t.every == 0) {
            return Err(CellweaveError::Config(
                "tracing every 0 steps: the clones are written down every so many steps, at \
                 least one"
                    .to_owned(),
            ));
        }
        if let Some(max) = spec.mechanics.and_then(|m| m.max_target)
            && (max.is_nan() || max < DIVISION_RATIO)
        {
            return Err(CellweaveError::Config(format!(
                "the network caps a cell's target at {max} times its reference size, and a \
                 cell divides at {DIVISION_RATIO}; no cell could ever divide"
            )));
        }
        let volume = scenario.depth > 1;
        if volume {
            refuse_in_a_volume(&scenario)?;
        }
        let mut lattice = CpmLattice::with_depth(
            scenario.width,
            scenario.height,
            scenario.depth,
            scenario.temperature,
        )
        .map_err(|e| CellweaveError::Config(e.to_string()))?;
        lattice.set_contact_order(scenario.contact_order);
        let cx = scenario.width / 2;
        let cy = scenario.height / 2;
        let cz = scenario.depth / 2;
        let first_radius = if volume {
            FIRST_CELL_RADIUS_3D
        } else {
            FIRST_CELL_RADIUS
        };

        match scenario.oxygen_source {
            OxygenSource::Bath => {}
            OxygenSource::Vessel { radius } => {
                lattice.add_obstacle_disc(Pos2 { x: cx, y: cy }, radius);
            }
            OxygenSource::Capillary { radius } => {
                for x in 0..scenario.width {
                    lattice.add_obstacle_disc(
                        Pos2 {
                            x,
                            y: capillary_row(scenario.height),
                        },
                        radius,
                    );
                }
            }
        }
        let first = match scenario.seeding {
            Seeding::Single => match scenario.oxygen_source {
                OxygenSource::Bath => lattice.add_ball(
                    CellKind::Tumor,
                    Pos3 {
                        x: cx,
                        y: cy,
                        z: cz,
                    },
                    first_radius,
                ),
                OxygenSource::Vessel { radius } => lattice.add_ball(
                    CellKind::Tumor,
                    Pos3 {
                        x: cx + radius + first_radius,
                        y: cy,
                        z: cz,
                    },
                    first_radius,
                ),
                OxygenSource::Capillary { radius } => lattice.add_cell(
                    CellKind::Tumor,
                    Pos2 {
                        x: cx,
                        y: capillary_row(scenario.height) + radius + FIRST_CELL_RADIUS,
                    },
                    FIRST_CELL_RADIUS,
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
                if scenario.oxygen_source.is_vessel() {
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
        let lesion_network = match &scenario.lesion {
            None => None,
            Some(lesion) => {
                let inside = |(x, y): (f64, f64)| {
                    (0.0..=f64::from(scenario.width)).contains(&x)
                        && (0.0..=f64::from(scenario.height)).contains(&y)
                };
                if lesion.centre.is_some_and(|c| !inside(c)) {
                    return Err(CellweaveError::Config(format!(
                        "a lesion centred at {:?} lies outside the {}x{} box",
                        lesion.centre, scenario.width, scenario.height
                    )));
                }
                if !(lesion.radius.is_finite() && lesion.radius > 0.0) {
                    return Err(CellweaveError::Config(format!(
                        "a lesion of radius {} holds no cell",
                        lesion.radius
                    )));
                }
                Some(Network::from_spec(&lesion.spec, newborn_volume)?)
            }
        };

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
            OxygenSource::Vessel { .. } | OxygenSource::Capillary { .. } => 0.0,
        };
        let oxygen = ScalarField::with_depth(
            scenario.width,
            scenario.height,
            scenario.depth,
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
                Some(ScalarField::with_depth(
                    scenario.width,
                    scenario.height,
                    scenario.depth,
                    STRESS_DIFFUSION,
                    0.0,
                    0.0,
                )?)
            }
        };

        let waste = match scenario.waste {
            None => None,
            Some(_) => Some(ScalarField::with_depth(
                scenario.width,
                scenario.height,
                scenario.depth,
                WASTE_DIFFUSION,
                0.0,
                0.0,
            )?),
        };
        let necrotic_factor = match scenario.necrotic_factor {
            None => None,
            Some(_) => Some(ScalarField::with_depth(
                scenario.width,
                scenario.height,
                scenario.depth,
                NECROTIC_FACTOR_DIFFUSION,
                0.0,
                0.0,
            )?),
        };
        let vegf = match scenario.angiogenesis {
            None => None,
            Some(a) => {
                if !scenario.oxygen_source.is_vessel() || scenario.oxygen_uptake == 0.0 {
                    return Err(CellweaveError::Config(
                        "vessels grow from a vessel, and this scenario has none that feeds \
                         anything; give it a vessel and oxygen uptake"
                            .into(),
                    ));
                }
                if !(a.threshold.is_finite() && a.threshold > 0.0) {
                    return Err(CellweaveError::Config(format!(
                        "the angiogenesis threshold is {}; it must be above zero, or every \
                         vessel sprouts whether anything calls or not",
                        a.threshold
                    )));
                }
                if !(a.reach.is_finite() && a.reach > 0.0) {
                    return Err(CellweaveError::Config(format!(
                        "VEGF's reach is {} pixels; it must be above zero",
                        a.reach
                    )));
                }
                if !(a.speed > 0.0 && a.speed <= 1.0) {
                    return Err(CellweaveError::Config(format!(
                        "a sprout's speed is {} pixels per step; it must be above zero and at \
                         most one, or the vessel it lays has gaps",
                        a.speed
                    )));
                }
                if a.regress_after == Some(0) {
                    return Err(CellweaveError::Config(
                        "a sprout that regresses after 0 steps is taken back the step it \
                         starts; give it at least one"
                            .into(),
                    ));
                }
                if a.perfusion == Perfusion::Flow
                    && !matches!(scenario.oxygen_source, OxygenSource::Capillary { .. })
                {
                    return Err(CellweaveError::Config(
                        "blood flows as pressure drives it from where it enters a capillary \
                         to where it leaves; give the scenario a capillary"
                            .into(),
                    ));
                }
                if let Some(floor) = a.prune_below {
                    if a.perfusion != Perfusion::Flow || a.regress_after.is_none() {
                        return Err(CellweaveError::Config(
                            "a sprout is pruned for carrying too little blood only when blood \
                             flows as pressure drives it and sprouts can regress; set \
                             perfusion to flow and regress_after"
                                .into(),
                        ));
                    }
                    if !(floor.is_finite() && floor > 0.0) {
                        return Err(CellweaveError::Config(format!(
                            "prune_below is {floor}; it must be above zero"
                        )));
                    }
                }
                if let Some(half) = a.half_delivery_flow {
                    if a.perfusion != Perfusion::Flow {
                        return Err(CellweaveError::Config(
                            "a vessel delivers in proportion to the blood through it only \
                             when blood flows as pressure drives it; set perfusion to flow"
                                .into(),
                        ));
                    }
                    if !(half.is_finite() && half > 0.0) {
                        return Err(CellweaveError::Config(format!(
                            "half_delivery_flow is {half}; it must be above zero"
                        )));
                    }
                }
                if a.regress_after.is_some() && a.perfusion == Perfusion::AtOnce {
                    return Err(CellweaveError::Config(
                        "a sprout regresses when it carries no blood, and with perfusion at \
                         once every sprout carries blood from its first pixel; set perfusion \
                         to on_joining"
                            .into(),
                    ));
                }
                if seed_network.angiogenic().is_none() {
                    return Err(CellweaveError::Config(
                        "vessels grow towards what cells secrete, and this network names no \
                         angiogenic node"
                            .into(),
                    ));
                }
                Some(ScalarField::new(
                    scenario.width,
                    scenario.height,
                    VEGF_DIFFUSION,
                    vegf_decay(&a),
                    0.0,
                )?)
            }
        };

        let chemokine = match scenario.immunity {
            None => None,
            Some(im) => {
                if !scenario.oxygen_source.is_vessel() || scenario.oxygen_uptake == 0.0 {
                    return Err(CellweaveError::Config(
                        "immune cells arrive through vessels carrying blood, and this \
                         scenario has none that feeds anything; give it a vessel and oxygen \
                         uptake"
                            .into(),
                    ));
                }
                let chance = |p: f64| p.is_finite() && p > 0.0 && p <= 1.0;
                if !chance(im.arrival) || !chance(im.kill) {
                    return Err(CellweaveError::Config(format!(
                        "immune arrival {} and kill {} are chances per step; each must be \
                         above zero and at most one",
                        im.arrival, im.kill
                    )));
                }
                if im.kills == 0 || im.lifetime == 0 {
                    return Err(CellweaveError::Config(
                        "an immune cell that can kill nothing, or lives no step, does nothing; \
                         give it kills and a lifetime"
                            .into(),
                    ));
                }
                if !(im.reach.is_finite() && im.reach > 0.0)
                    || !(im.chemotaxis.is_finite() && im.chemotaxis >= 0.0)
                {
                    return Err(CellweaveError::Config(format!(
                        "the chemokine's reach is {} pixels and the pull up it {}; the reach \
                         must be above zero and the pull not below it",
                        im.reach, im.chemotaxis
                    )));
                }
                Some(ScalarField::new(
                    scenario.width,
                    scenario.height,
                    VEGF_DIFFUSION,
                    VEGF_DIFFUSION / (im.reach * im.reach),
                    0.0,
                )?)
            }
        };

        let mut vessels = Vessels::default();
        if let OxygenSource::Capillary { radius } = scenario.oxygen_source {
            let flowing = scenario
                .angiogenesis
                .is_some_and(|a| a.perfusion == Perfusion::Flow);
            vessels.lie_along(
                Capillary {
                    row: capillary_row(scenario.height),
                    width: scenario.width,
                    radius,
                },
                flowing,
            );
        }

        Ok(Self {
            rng: StdRng::seed_from_u64(scenario.seed),
            scenario,
            lattice,
            networks,
            lineage,
            clones,
            pristine: seed_network,
            lives: Lives::new(cells),
            lesion_network,
            oxygen,
            stress,
            stress_level: 0.0,
            vegf,
            waste,
            necrotic_factor,
            vessels,
            first,
            next_mcs: 0,
            born: cells as u64,
            gone: 0,
            crushed: 0,
            chemokine,
            immune: Vec::new(),
            killed: 0,
            marks: vec![None; cells],
            kept_dead: BTreeSet::new(),
            died_at: BTreeMap::new(),
            debris_owed: 0.0,
            bath_level: 1.0,
            due: vec![0; cells],
            began: vec![0; cells],
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
        self.settle_bath_level();
        let running = Simulation {
            scenario: self.scenario,
            lattice: self.lattice,
            networks: self.networks,
            lineage: self.lineage,
            clones: self.clones,
            pristine: self.pristine,
            lives: self.lives,
            lesion_network: self.lesion_network,
            oxygen: self.oxygen,
            stress: self.stress,
            stress_level: self.stress_level,
            vegf: self.vegf,
            waste: self.waste,
            necrotic_factor: self.necrotic_factor,
            vessels: self.vessels,
            rng: self.rng,
            first: self.first,
            next_mcs: self.next_mcs,
            born: self.born,
            gone: self.gone,
            crushed: self.crushed,
            chemokine: self.chemokine,
            immune: self.immune,
            killed: self.killed,
            marks: self.marks,
            kept_dead: self.kept_dead,
            died_at: self.died_at,
            debris_owed: self.debris_owed,
            bath_level: self.bath_level,
            due: self.due,
            began: self.began,
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
        if self.scenario.revert.as_ref().is_some_and(|r| r.at == mcs) {
            self.revert_clones(mcs);
        }
        if self.scenario.lesion.as_ref().is_some_and(|l| l.at == mcs) {
            self.place_lesion(mcs);
        }
        if self.scenario.tracing.is_some_and(|t| t.at == mcs) {
            self.mark_every_cell();
        }

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
        self.settle_bath_level();
        let local_o2 = mean_per_cell(&self.oxygen, &self.lattice);
        let exposure = self.settle_stress(mcs)?;
        let local_stress = self
            .stress
            .as_ref()
            .map(|field| mean_per_cell(field, &self.lattice));
        self.settle_waste(mcs)?;
        let local_waste = self
            .waste
            .as_ref()
            .map(|field| mean_per_cell(field, &self.lattice));
        self.settle_necrotic_factor(mcs)?;
        let local_factor = self
            .necrotic_factor
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
                waste: Concentration(local_waste.as_ref().map_or(0.0, |w| w[i + 1])),
                necrotic_factor: Concentration(local_factor.as_ref().map_or(0.0, |f| f[i + 1])),
            };
            for _ in 0..NETWORK_TIME_PER_MCS {
                net.step(TimeStep(1.0), &inputs);
            }
        }

        // 2b. With vessels that can grow, what the living secrete settles into
        //     the VEGF field, and the vessels answer it. Laying vessel takes
        //     pixels from cells, and a cell left with none was crushed.
        self.call_for_vessels(mcs)?;

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
                    if volume > 0 && leaving == CellKind::Necrotic {
                        self.died_at.entry(i).or_insert(mcs);
                    }
                    if volume == 0 {
                        *slot = None;
                        self.gone += 1;
                        self.kept_dead.remove(&i);
                        self.died_at.remove(&i);
                    } else if let (CellKind::Necrotic, Some(loss)) =
                        (leaving, self.scenario.debris_loss)
                    {
                        // Its debris breaks down where it lies: the cell
                        // keeps to what it held when it died (below, or here
                        // if it was made necrotic outside this loop), and the
                        // dead tissue owes the fluid a share of that each
                        // step.
                        if self.kept_dead.insert(i) {
                            self.lattice.set_volume_target(
                                id,
                                Volume(volume),
                                net.mechanics().lambda_volume,
                            );
                        }
                        self.debris_owed += loss * f64::from(volume);
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
                    self.lives.left(i);
                    self.died_at.insert(i, mcs);
                    // With a debris loss its debris is the dead tissue's,
                    // and owes its share, from this step. Kept from the
                    // next step only, the cells dead one step were neither
                    // living nor kept, a shell between the kept core and
                    // the living, and nothing closed over the core one step
                    // in two.
                    if let Some(loss) = self.scenario.debris_loss
                        && self.kept_dead.insert(i)
                    {
                        let volume = self.lattice.volume(id).map_or(0, |v| v.0);
                        self.lattice.set_volume_target(
                            id,
                            Volume(volume),
                            net.mechanics().lambda_volume,
                        );
                        self.debris_owed += loss * f64::from(volume);
                    }
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
                    self.lives.left(i);
                }
                Some(_)
                    if net
                        .differentiation()
                        .is_some_and(|chance| self.rng.random::<f64>() < chance) =>
                {
                    self.lattice.set_kind(id, CellKind::Differentiated);
                    self.lives.left(i);
                }
                Some(_) => {
                    let m = net.mechanics();
                    // With a cycle a cell grows through it, rather than
                    // asking at birth for all it will be.
                    let target = match self.scenario.cycle {
                        Some(_) => grown(m.target_volume, self.began[i], self.due[i], mcs),
                        None => m.target_volume,
                    };
                    self.lattice.set_volume_target(id, target, m.lambda_volume);
                }
                None => {}
            }
        }

        // 3a. With debris loss, the dead tissue gives up at its edge every
        //     whole pixel it owes the fluid, and fluid shut in the tissue
        //     drains where it touches the living.
        self.break_down_debris();
        self.drain_shut_in();

        // 3b. Immune cells arrive, kill what they touch, and die; then the
        //     chemokine settles and they are set to climb it.
        self.immune_respond(mcs)?;

        // 4. Advance CPM.
        for _ in 0..self.scenario.sweeps_per_step {
            self.lattice.monte_carlo_step(&mut self.rng);
        }

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
        let every = match (self.scenario.seeding, self.scenario.cycle) {
            (Seeding::Single, None) => DIVISION_CHECK_EVERY,
            _ => 1,
        };
        if mcs.is_multiple_of(every) && self.living() < self.scenario.max_cells {
            self.born += self.divide_ready_cells(mcs);
        }

        Ok(report)
    }

    /// Settle the VEGF field against what every living cell secretes, and
    /// let the vessels grow. Nothing happens without angiogenesis.
    fn call_for_vessels(&mut self, mcs: u64) -> Result<(), CellweaveError> {
        let (Some(angiogenesis), Some(vegf)) = (self.scenario.angiogenesis, self.vegf.as_mut())
        else {
            return Ok(());
        };
        // Secreted per pixel at the decay rate times the node's level, so
        // that tissue all secreting one level settles at that level.
        let decay = vegf_decay(&angiogenesis);
        let mut per_cell = vec![0.0; self.networks.len() + 1];
        for (i, slot) in self.networks.iter().enumerate() {
            let Some(net) = slot else { continue };
            let alive = matches!(
                self.lattice.cell_kind(cell_of(i)),
                Some(CellKind::Tumor | CellKind::Normal)
            );
            if alive {
                per_cell[i + 1] = decay * net.angiogenic().unwrap_or(0.0);
            }
        }
        apply_production(vegf, &self.lattice, &per_cell);
        vegf.relax(TimeStep(1.0), VEGF_TOLERANCE * decay, OXYGEN_MAX_STEPS)
            .map_err(|e| CellweaveError::Simulation(format!("VEGF at MCS {mcs}: {e}")))?;
        if self
            .vessels
            .grow(&mut self.lattice, vegf, &angiogenesis, mcs)
        {
            for i in 0..self.networks.len() {
                let id = cell_of(i);
                let crushed = self.networks[i].is_some()
                    && matches!(
                        self.lattice.cell_kind(id),
                        Some(CellKind::Tumor | CellKind::Normal)
                    )
                    && self.lattice.volume(id).is_some_and(|v| v.0 == 0);
                if crushed {
                    self.lattice.set_kind(id, CellKind::Necrotic);
                    self.lives.left(i);
                    self.crushed += 1;
                }
            }
        }
        Ok(())
    }

    /// With immunity: maybe one immune cell arrives through a vessel wall
    /// carrying blood; each immune cell touching a living tumour cell may
    /// kill one; one exhausted or past its lifetime dies and breaks up; and
    /// the chemokine living tumour cells secrete settles, for the immune
    /// cells to climb. Nothing happens, and nothing is drawn, without it.
    fn immune_respond(&mut self, mcs: u64) -> Result<(), CellweaveError> {
        let (Some(immunity), Some(_)) = (self.scenario.immunity, self.chemokine.as_ref()) else {
            return Ok(());
        };
        if mcs < immunity.from {
            return Ok(());
        }
        if self.rng.random::<f64>() < immunity.arrival {
            self.immune_arrives(immunity, mcs);
        }
        self.immune_kill(immunity, mcs);
        let decay = VEGF_DIFFUSION / (immunity.reach * immunity.reach);
        let mut per_cell = vec![0.0; self.networks.len() + 1];
        for (i, slot) in self.networks.iter().enumerate() {
            if slot.is_some() && self.lattice.cell_kind(cell_of(i)) == Some(CellKind::Tumor) {
                per_cell[i + 1] = decay;
            }
        }
        let Some(field) = self.chemokine.as_mut() else {
            return Ok(());
        };
        apply_production(field, &self.lattice, &per_cell);
        field
            .relax(TimeStep(1.0), VEGF_TOLERANCE * decay, OXYGEN_MAX_STEPS)
            .map_err(|e| CellweaveError::Simulation(format!("chemokine at MCS {mcs}: {e}")))?;
        let (w, h) = self.lattice.dims();
        let values = (0..h)
            .flat_map(|y| (0..w).map(move |x| (x, y)))
            .map(|(x, y)| field.concentration_at(x, y))
            .collect();
        self.lattice.attract_immune(values, immunity.chemotaxis);
        Ok(())
    }

    /// One immune cell arrives beside a vessel wall carrying blood, at a
    /// place drawn among all such places; none if there is none.
    fn immune_arrives(&mut self, immunity: Immunity, mcs: u64) {
        let lattice = &self.lattice;
        let (w, h) = lattice.dims();
        let beside_blood = |x: u32, y: u32| {
            lattice.occupant(x, y) != OBSTACLE
                && [
                    (x.wrapping_sub(1), y),
                    (x + 1, y),
                    (x, y.wrapping_sub(1)),
                    (x, y + 1),
                ]
                .iter()
                .any(|&(i, j)| {
                    i < w
                        && j < h
                        && lattice.occupant(i, j) == OBSTACLE
                        && self.vessels.carries_blood(i, j)
                })
        };
        let places: Vec<(u32, u32)> = (0..h)
            .flat_map(|y| (0..w).map(move |x| (x, y)))
            .filter(|&(x, y)| beside_blood(x, y))
            .collect();
        if places.is_empty() {
            return;
        }
        let pick = (self.rng.random::<f64>() * places.len() as f64) as usize;
        let (x, y) = places[pick.min(places.len() - 1)];
        let Some(id) = self
            .lattice
            .seed_cell(CellKind::Immune, Pos2 { x, y }, IMMUNE_RADIUS)
        else {
            return;
        };
        self.lattice
            .set_volume_target(id, Volume(IMMUNE_VOLUME), IMMUNE_LAMBDA);
        self.networks.push(None);
        self.lineage.push(CloneId(0));
        self.marks.push(None);
        self.due.push(0);
        self.began.push(0);
        self.immune.push(ImmuneCell {
            slot: self.networks.len() - 1,
            arrived: mcs,
            kills_left: immunity.kills,
            dying: false,
        });
        // A cell it arrived on top of may have been left with nothing; it
        // dies of it, though no vessel crushed it.
        for i in 0..self.networks.len() {
            let id = cell_of(i);
            let emptied = self.networks[i].is_some()
                && matches!(
                    self.lattice.cell_kind(id),
                    Some(CellKind::Tumor | CellKind::Normal)
                )
                && self.lattice.volume(id).is_some_and(|v| v.0 == 0);
            if emptied {
                self.lattice.set_kind(id, CellKind::Necrotic);
                self.lives.left(i);
            }
        }
    }

    /// Each immune cell touching a living tumour cell may kill the one of
    /// them with the lowest identifier not killed already this step: it
    /// starts apoptosis. One exhausted or past its lifetime dies and breaks
    /// up, and one broken up, or squeezed out of every pixel, is gone.
    fn immune_kill(&mut self, immunity: Immunity, mcs: u64) {
        if self.immune.is_empty() {
            return;
        }
        let lattice = &self.lattice;
        let (w, h) = lattice.dims();
        let mut which = vec![usize::MAX; lattice.cell_count() + 1];
        for (k, cell) in self.immune.iter().enumerate() {
            if !cell.dying {
                which[cell.slot + 1] = k;
            }
        }
        let mut touching: Vec<Vec<u32>> = vec![Vec::new(); self.immune.len()];
        for y in 0..h {
            for x in 0..w {
                let id = lattice.occupant(x, y);
                let Some(&k) = which.get(id as usize).filter(|&&k| k != usize::MAX) else {
                    continue;
                };
                for (i, j) in [
                    (x + 1, y),
                    (x, y + 1),
                    (x.wrapping_sub(1), y),
                    (x, y.wrapping_sub(1)),
                ] {
                    if i >= w || j >= h {
                        continue;
                    }
                    let other = lattice.occupant(i, j);
                    let alive = matches!(
                        self.fate(CellId(other)),
                        Some(Fate::Cycling | Fate::Quiescent)
                    ) && lattice.cell_kind(CellId(other)) == Some(CellKind::Tumor);
                    if other != id
                        && other != 0
                        && other != OBSTACLE
                        && alive
                        && !touching[k].contains(&other)
                    {
                        touching[k].push(other);
                    }
                }
            }
        }
        for (k, targets) in touching.iter_mut().enumerate() {
            targets.sort_unstable();
            let cell = self.immune[k];
            let id = cell_of(cell.slot);
            if cell.dying {
                let volume = self.lattice.volume(id).map_or(0, |v| v.0);
                let target = Volume(volume.saturating_sub(FRAGMENTATION_RATE));
                self.lattice.set_volume_target(id, target, IMMUNE_LAMBDA);
                continue;
            }
            if cell.kills_left == 0 || mcs.saturating_sub(cell.arrived) >= immunity.lifetime {
                self.immune[k].dying = true;
                self.lattice.set_kind(id, CellKind::Apoptotic);
                continue;
            }
            // Another immune cell may have killed it earlier this step.
            let Some(&victim) = targets
                .iter()
                .find(|&&v| self.lattice.cell_kind(CellId(v)) == Some(CellKind::Tumor))
            else {
                continue;
            };
            // A tumour cell whose network names an evasion node is shielded
            // by its level.
            let shield = self.networks[victim as usize - 1]
                .as_ref()
                .and_then(Network::evasion)
                .map_or(0.0, |e| e.clamp(0.0, 1.0));
            if self.rng.random::<f64>() < immunity.kill * (1.0 - shield) {
                self.lattice.set_kind(CellId(victim), CellKind::Apoptotic);
                self.lives.left(victim as usize - 1);
                self.immune[k].kills_left -= 1;
                self.killed += 1;
            }
        }
        // One with no pixel left, broken up or squeezed out, is gone.
        let lattice = &self.lattice;
        self.immune
            .retain(|c| lattice.volume(cell_of(c.slot)).is_some_and(|v| v.0 > 0));
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
                big_enough && self.due[*i] <= mcs && self.fate(id) == Some(Fate::Cycling)
            })
            .map(|(i, _)| i)
            .collect();
        if ready.is_empty() {
            return 0;
        }

        let mut born = 0;
        let living = self.living();
        let ids: Vec<CellId> = ready.iter().map(|&i| cell_of(i)).collect();
        let pixels = self.lattice.pixels_of(&ids);
        for (i, pixels) in ready.into_iter().zip(pixels) {
            if living + born as usize >= self.scenario.max_cells {
                break;
            }
            // In the plane, a line at an angle drawn uniformly; in a volume, a
            // plane whose normal is drawn uniformly over the sphere, with two
            // draws, so that a run in the plane draws what it always did.
            let normal = if self.scenario.depth > 1 {
                let z = 2.0 * self.rng.random::<f64>() - 1.0;
                let around = self.rng.random::<f64>() * std::f64::consts::TAU;
                let r = (1.0 - z * z).max(0.0).sqrt();
                (r * around.cos(), r * around.sin(), z)
            } else {
                let angle = self.rng.random::<f64>() * std::f64::consts::PI;
                (angle.cos(), angle.sin(), 0.0)
            };
            if self
                .lattice
                .divide_across(cell_of(i), &pixels, normal)
                .is_some()
            {
                let mut daughter = self.networks[i].clone();
                let mut clone = self.lineage[i];
                // The rate is tested before anything is drawn, so that a run
                // without mutation takes exactly the draws it always took.
                let rate = self.scenario.mutation_rate;
                if rate > 0.0 && self.rng.random::<f64>() < rate {
                    let (which, how) = (self.rng.random::<f64>(), self.rng.random::<f64>());
                    // Drawn only when rewiring is on, for the same reason as
                    // the mutation rate: no rewiring, no extra draw.
                    let rewire = self.scenario.rewiring;
                    let wiring = if rewire > 0.0 {
                        let u = self.rng.random::<f64>();
                        (u < rewire).then(|| {
                            let (source, sign) =
                                (self.rng.random::<f64>(), self.rng.random::<f64>());
                            [u / rewire, which, source, sign, how]
                        })
                    } else {
                        None
                    };
                    let change = |n: &mut Network| match wiring {
                        Some(draws) => n.rewire(draws),
                        None => n.mutate(which, how),
                    };
                    let founded = CloneId(self.clones.len() as u32);
                    let suppressed = self.scenario.suppress.contains(&founded)
                        || self
                            .scenario
                            .suppress_before
                            .is_some_and(|before| mcs < before);
                    // A suppressed mutation is made on a copy, to know what it
                    // would have been, and the copy is thrown away.
                    let made = if suppressed {
                        daughter.clone().and_then(|mut n| change(&mut n))
                    } else {
                        daughter.as_mut().and_then(change)
                    };
                    if let Some(mutation) = made {
                        self.clones.push(CloneRecord {
                            id: founded,
                            parent: Some(clone),
                            founded_at: mcs,
                            mutation: Some(mutation),
                            suppressed,
                            lesion: false,
                            reverted_at: None,
                        });
                        clone = founded;
                    }
                }
                // Drawn last, and only with a cycle, so that a run without
                // one takes exactly the draws it always took.
                let (began, due) = match self.scenario.cycle {
                    Some(cycle) => {
                        self.due[i] = mcs + self.cycle_length(cycle);
                        self.began[i] = mcs;
                        (mcs, mcs + self.cycle_length(cycle))
                    }
                    None => (0, 0),
                };
                self.networks.push(daughter);
                self.lineage.push(clone);
                self.marks.push(self.marks[i]);
                self.due.push(due);
                self.began.push(began);
                self.lives.divided(i, self.networks.len() - 1);
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

    /// How many steps a cycle just begun lasts, drawn from the run's
    /// generator: the median times $e^{\sigma z}$, with $z$ a standard normal
    /// draw (Box and Muller, two uniform draws) and $\sigma$ the spread. One
    /// step at least.
    fn cycle_length(&mut self, cycle: Cycle) -> u64 {
        let (u, v) = (1.0 - self.rng.random::<f64>(), self.rng.random::<f64>());
        let z = (-2.0 * u.ln()).sqrt() * (std::f64::consts::TAU * v).cos();
        (cycle.median * (cycle.spread * z).exp()).round().max(1.0) as u64
    }

    /// Take every whole pixel the dead tissue's debris owes off it at its
    /// edge, and close the tissue over each (`CpmLattice::close_over`): the
    /// pixels are drawn among dead ones that touch living tissue, and once
    /// those are taken, the dead pixels they laid bare, up to
    /// `EDGE_PASSES` times. A dead cell inside the core has only dead cells
    /// around it, and a pixel it gave one of them would only move the debris
    /// about; at the edge the living tissue closes in, which is what the
    /// core's loss comes to wherever in it the debris breaks down. What is
    /// still owed is carried to the next step, never more than the dead
    /// tissue holds. Draws nothing without debris loss.
    fn break_down_debris(&mut self) {
        if self.scenario.debris_loss.is_none() {
            return;
        }
        let owed = self.debris_owed.floor() as usize;
        if owed == 0 {
            return;
        }
        let slots: Vec<usize> = self.kept_dead.iter().copied().collect();
        let ids: Vec<CellId> = slots.iter().map(|&i| cell_of(i)).collect();
        let mut pool: Vec<(usize, usize)> = self
            .lattice
            .pixels_of(&ids)
            .into_iter()
            .enumerate()
            .flat_map(|(k, pixels)| pixels.into_iter().map(move |p| (k, p)))
            .collect();
        let held = pool.len();
        let mut given = vec![0u32; slots.len()];
        let mut paid = 0;
        for _ in 0..EDGE_PASSES {
            let (mut edge, inside): (Vec<_>, Vec<_>) = pool
                .into_iter()
                .partition(|&(_, p)| self.lattice.can_close_over(p));
            pool = inside;
            if edge.is_empty() {
                break;
            }
            while paid < owed && !edge.is_empty() {
                let (k, p) = edge.swap_remove(self.rng.random_range(0..edge.len()));
                if self.lattice.close_over(p, &mut self.rng) {
                    given[k] += 1;
                    paid += 1;
                }
            }
            if paid == owed {
                break;
            }
            pool.extend(edge);
        }
        self.debris_owed = (self.debris_owed - paid as f64).min((held - paid) as f64);
        // What is lost is lost: each dead cell tries to hold that much less,
        // and does not grow back into the space it gave up.
        for (&i, lost) in slots.iter().zip(given) {
            let id = cell_of(i);
            if let (true, Some(target), Some(net)) =
                (lost > 0, self.lattice.target_volume(id), &self.networks[i])
            {
                let target = Volume(target.0.saturating_sub(lost));
                self.lattice
                    .set_volume_target(id, target, net.mechanics().lambda_volume);
            }
        }
    }

    /// With a debris loss, close in every medium pixel the tissue shuts in
    /// that touches living tissue, one layer a step: fluid in a cavity of
    /// the tissue leaves through it as the debris does, and the living cells
    /// around it take the space (`CpmLattice::close_in`). Which neighbour
    /// takes each is drawn from the run's generator, in the order of the
    /// pixels. Without a debris loss nothing is looked at or drawn.
    fn drain_shut_in(&mut self) {
        if self.scenario.debris_loss.is_none() {
            return;
        }
        let edge: Vec<usize> = self
            .lattice
            .shut_in()
            .into_iter()
            .filter(|&i| self.lattice.can_close_in(i))
            .collect();
        for i in edge {
            self.lattice.close_in(i, &mut self.rng);
        }
    }

    /// Undo the founding mutation of every clone the reversion names, in every
    /// cell whose lineage passes through one: its network is rebuilt from the
    /// network its lineage started from, with every mutation made since but
    /// those, and keeps the protein levels it had. Draws nothing.
    fn revert_clones(&mut self, mcs: u64) {
        let Some(undone) = self.scenario.revert.as_ref().map(|r| r.clones.clone()) else {
            return;
        };
        for id in &undone {
            if let Some(record) = self.clones.get_mut(id.0 as usize)
                && record.mutation.is_some()
            {
                record.reverted_at = Some(mcs);
            }
        }
        let mut touched: Vec<bool> = Vec::with_capacity(self.clones.len());
        for clone in &self.clones {
            let inherited = clone
                .parent
                .and_then(|p| touched.get(p.0 as usize).copied())
                .unwrap_or(false);
            touched.push(inherited || undone.contains(&clone.id));
        }
        let rebuilt: Vec<(usize, Network)> = self
            .networks
            .iter()
            .enumerate()
            .filter(|&(i, _)| {
                touched
                    .get(self.lineage[i].0 as usize)
                    .copied()
                    .unwrap_or(false)
            })
            .filter_map(|(i, slot)| {
                let current = slot.as_ref()?;
                let mut network = self.lineage_network(self.lineage[i], &undone);
                network.take_state(current);
                Some((i, network))
            })
            .collect();
        for (i, network) in rebuilt {
            self.networks[i] = Some(network);
        }
    }

    /// The network of a clone's lineage, rebuilt: the network it started
    /// from, the unmutated one or a lesion's, with every mutation made since,
    /// in order, but those of the clones in `leaving_out`.
    fn lineage_network(&self, clone: CloneId, leaving_out: &[CloneId]) -> Network {
        let mut chain = Vec::new();
        let mut at = self.clones.get(clone.0 as usize);
        while let Some(record) = at {
            chain.push(record);
            if record.lesion {
                break;
            }
            at = record.parent.and_then(|p| self.clones.get(p.0 as usize));
        }
        chain.reverse();
        let base = match (chain.first(), &self.lesion_network) {
            (Some(first), Some(lesion)) if first.lesion => lesion,
            _ => &self.pristine,
        };
        let mut network = base.clone();
        for record in chain {
            if record.suppressed || record.reverted_at.is_some() || leaving_out.contains(&record.id)
            {
                continue;
            }
            if let Some(mutation) = &record.mutation {
                // A mutation that no longer fits, because one it relied on
                // was left out, is left out too.
                let _ = network.apply(mutation);
            }
        }
        network
    }

    /// Mark every cell in the layer, each with its own slot, as the basal
    /// cells an experiment marks are of every kind. Draws nothing.
    fn mark_every_cell(&mut self) {
        for i in 0..self.marks.len() {
            let in_layer = self.in_layer(cell_of(i));
            if let Some(mark) = self.marks.get_mut(i) {
                *mark = in_layer.then_some(i as u32);
            }
        }
    }

    /// Turn every living cell whose centre lies within the lesion's radius of
    /// the centre of the box into a lesion cell: the lesion's network, fresh,
    /// and one new clone for all of them. Draws nothing.
    fn place_lesion(&mut self, mcs: u64) {
        let (Some(lesion), Some(network)) = (&self.scenario.lesion, &self.lesion_network) else {
            return;
        };
        let (w, h) = self.lattice.dims();
        let (cx, cy) = lesion
            .centre
            .unwrap_or((f64::from(w) / 2.0, f64::from(h) / 2.0));
        // Centres of mass, one pass over the lattice, indexed by cell.
        let mut sums = vec![(0.0_f64, 0.0_f64, 0_u32); self.lattice.cell_count() + 1];
        for y in 0..h {
            for x in 0..w {
                let id = self.lattice.occupant(x, y);
                if let Some(s) = sums.get_mut(id as usize) {
                    s.0 += f64::from(x) + 0.5;
                    s.1 += f64::from(y) + 0.5;
                    s.2 += 1;
                }
            }
        }
        let joining: Vec<usize> = (0..self.networks.len())
            .filter(|&i| {
                let id = cell_of(i);
                let Some(&(sx, sy, n)) = sums.get(id.0 as usize) else {
                    return false;
                };
                let alive = matches!(self.fate(id), Some(Fate::Cycling | Fate::Quiescent));
                alive
                    && n > 0
                    && (sx / f64::from(n) - cx).hypot(sy / f64::from(n) - cy) <= lesion.radius
            })
            .collect();
        if joining.is_empty() {
            return;
        }
        let founded = CloneId(self.clones.len() as u32);
        self.clones.push(CloneRecord {
            id: founded,
            parent: None,
            founded_at: mcs,
            mutation: None,
            suppressed: false,
            lesion: true,
            reverted_at: None,
        });
        for i in joining {
            self.networks[i] = Some(network.clone());
            self.lineage[i] = founded;
        }
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
            OxygenSource::Vessel { .. } | OxygenSource::Capillary { .. } => {
                apply_background_uptake(field, &self.lattice, exposure.uptake);
                match self
                    .scenario
                    .angiogenesis
                    .and_then(|a| a.half_delivery_flow)
                {
                    None => hold_vessels_at(field, &self.lattice, level, |x, y| {
                        self.vessels.carries_blood(x, y)
                    }),
                    Some(half) => hold_vessels_by(field, &self.lattice, level, |x, y| {
                        delivered(&self.vessels, x, y, half)
                    }),
                }
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

    /// Solve the waste field for step `mcs`: every pixel of a living cell
    /// releases the scenario's production, the dead and the medium release
    /// nothing, and the medium the bath reaches is held at zero. Without a
    /// waste there is no field and nothing is solved.
    fn settle_waste(&mut self, mcs: u64) -> Result<(), CellweaveError> {
        let (Some(waste), Some(field)) = (self.scenario.waste, self.waste.as_mut()) else {
            return Ok(());
        };
        settle_released(
            field,
            &self.lattice,
            self.networks.len(),
            waste.production,
            |_, kind| matches!(kind, CellKind::Tumor | CellKind::Normal),
        )
        .map_err(|e| CellweaveError::Simulation(format!("waste at MCS {mcs}: {e}")))
    }

    /// Solve the necrotic factor's field for step `mcs`: every pixel of a
    /// dead cell releases the scenario's production, for as many steps after
    /// the one it died at as the scenario says or for as long as it is there,
    /// the living and the medium release nothing, and the medium the bath
    /// reaches is held at zero. Without a necrotic factor there is no field
    /// and nothing is solved.
    fn settle_necrotic_factor(&mut self, mcs: u64) -> Result<(), CellweaveError> {
        let (Some(factor), Some(field)) =
            (self.scenario.necrotic_factor, self.necrotic_factor.as_mut())
        else {
            return Ok(());
        };
        let died_at = &self.died_at;
        let releasing = |id: CellId| match (factor.release_for, died_at.get(&(id.0 as usize - 1))) {
            (Some(steps), Some(&died)) => mcs <= died + steps,
            _ => true,
        };
        settle_released(
            field,
            &self.lattice,
            self.networks.len(),
            factor.production,
            |id, kind| kind == CellKind::Necrotic && releasing(id),
        )
        .map_err(|e| CellweaveError::Simulation(format!("necrotic factor at MCS {mcs}: {e}")))
    }

    /// With a still medium, the level the bath holds next to the tissue at
    /// the next step, from what the tissue took up at this one.
    ///
    /// A still medium held at the source's level far away and taking up
    /// nothing itself sits, around a ball of radius $R$ that takes up $F$ in
    /// all, at $1 - F / (4 \pi D R)$ at its surface. $F$ is what the held
    /// medium handed the tissue, read from the settled field; the ball is the
    /// one holding everything the bath does not reach, the tissue and the
    /// fluid shut in it. The level lags the field by a step: what a tissue
    /// takes up is set by its living cells, which a step changes by a few
    /// in a hundred, and a cell starved by a lower level takes many steps to
    /// die.
    fn settle_bath_level(&mut self) {
        let Some(medium) = self.scenario.still_medium else {
            return;
        };
        let taken = self.oxygen.inflow();
        let radius = (3.0 * self.oxygen.free_count() as f64 / (4.0 * std::f64::consts::PI)).cbrt();
        self.bath_level = if radius > 0.0 {
            (1.0 - taken / (4.0 * std::f64::consts::PI * medium.diffusion * radius)).max(0.0)
        } else {
            1.0
        };
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
            OxygenSource::Bath => {
                hold_medium_bath(&mut self.oxygen, &self.lattice, self.bath_level);
            }
            OxygenSource::Vessel { .. } | OxygenSource::Capillary { .. } => {
                apply_background_uptake(&mut self.oxygen, &self.lattice, uptake);
                match self
                    .scenario
                    .angiogenesis
                    .and_then(|a| a.half_delivery_flow)
                {
                    None => hold_vessels_at(&mut self.oxygen, &self.lattice, 1.0, |x, y| {
                        self.vessels.carries_blood(x, y)
                    }),
                    Some(half) => hold_vessels_by(&mut self.oxygen, &self.lattice, 1.0, |x, y| {
                        delivered(&self.vessels, x, y, half)
                    }),
                }
            }
        }
    }

    /// Whether any cell occupies a pixel on the edge of the box.
    fn touches_the_edge(&self) -> bool {
        let lattice = &self.lattice;
        let (w, h, d) = lattice.dims3();
        let slice = w as usize * h as usize;
        // Only a cell carrying a network is tissue: an immune cell wandering
        // to the edge, alive or breaking up, is not the tissue reaching it.
        let is_cell = |x: u32, y: u32, z: u32| {
            let id = lattice.occupant_at(z as usize * slice + (y * w + x) as usize);
            id != 0
                && id != OBSTACLE
                && self
                    .networks
                    .get(id as usize - 1)
                    .is_some_and(Option::is_some)
        };
        // Along a capillary the left and right sides of the box cut through a
        // layer that goes on beyond them, the way the plane cuts through a
        // cord; only the edges away from the capillary bound the tissue. In a
        // volume, a vessel runs through the front and back faces, which cut
        // through a cord the same way.
        let sides = !matches!(self.scenario.oxygen_source, OxygenSource::Capillary { .. });
        let faces = d > 1 && !self.scenario.oxygen_source.is_vessel();
        (0..d).any(|z| {
            (0..w).any(|x| is_cell(x, 0, z) || is_cell(x, h - 1, z))
                || (sides && (0..h).any(|y| is_cell(0, y, z) || is_cell(w - 1, y, z)))
        }) || (faces && (0..h).any(|y| (0..w).any(|x| is_cell(x, y, 0) || is_cell(x, y, d - 1))))
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

    /// The waste field, as it was last settled, or `None` when the scenario
    /// has no waste.
    pub fn waste(&self) -> Option<&ScalarField> {
        self.waste.as_ref()
    }

    /// The necrotic factor's field, as it was last settled, or `None` when
    /// the scenario has none.
    pub fn necrotic_factor(&self) -> Option<&ScalarField> {
        self.necrotic_factor.as_ref()
    }

    /// The level the bath holds the medium at next to the tissue, this step:
    /// 1, the source's, unless the medium is still.
    #[must_use]
    pub fn bath_level(&self) -> f64 {
        self.bath_level
    }

    /// What the living secrete as a call for vessels, when vessels can grow.
    pub fn vegf(&self) -> Option<&ScalarField> {
        self.vegf.as_ref()
    }

    /// The tips of the sprouts still growing.
    pub fn sprouts(&self) -> &[Tip] {
        self.vessels.tips()
    }

    /// How many sprouts have started, growing or not, standing or not.
    pub fn sprouted(&self) -> usize {
        self.vessels.sprouted()
    }

    /// How many sprouts have regressed, for carrying no blood too long.
    pub fn regressed(&self) -> usize {
        self.vessels.regressed()
    }

    /// The flow through the vessel at `(x, y)`, as a share of what the
    /// capillary alone carries, when blood flows as pressure drives it; zero
    /// otherwise, and meaningless for a pixel that is not vessel.
    pub fn flow_through(&self, x: u32, y: u32) -> f64 {
        self.vessels.flow_through(x, y)
    }

    /// Whether blood flows as pressure drives it, so that a vessel carries
    /// some share of it rather than all or none.
    pub fn flowing(&self) -> bool {
        self.scenario
            .angiogenesis
            .is_some_and(|a| a.perfusion == Perfusion::Flow)
    }

    /// Whether the vessel at `(x, y)` carries blood, and with it oxygen:
    /// the vessel the scenario laid always does, a sprout once perfused.
    /// Meaningless for a pixel that is not vessel.
    pub fn carries_blood(&self, x: u32, y: u32) -> bool {
        self.vessels.carries_blood(x, y)
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

    /// Pairs of sisters, by what both did next: divided again or left. Only
    /// pairs whose two next events have both happened are counted.
    pub fn sister_fates(&self) -> SisterFates {
        self.lives.fates()
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
                    CellKind::Necrotic
                        | CellKind::Differentiated
                        | CellKind::Apoptotic
                        | CellKind::Immune
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

    /// The immune cells in the tissue, not counting those dying.
    pub fn immune_cells(&self) -> usize {
        self.immune.iter().filter(|c| !c.dying).count()
    }

    /// Tumour cells immune cells have killed so far.
    pub fn killed(&self) -> u64 {
        self.killed
    }

    /// Whether a cell is in the layer as an experiment counts basal cells:
    /// cycling, quiescent, or committed to differentiate and not yet gone.
    fn in_layer(&self, id: CellId) -> bool {
        matches!(
            self.fate(id),
            Some(Fate::Cycling | Fate::Quiescent | Fate::Differentiating)
        ) && self.lattice.volume(id).is_some_and(|v| v.0 > 0)
    }

    /// The cells in the layer carrying each mark, as `(mark, cells)`, for
    /// every mark some cell in the layer still carries, in the order of the
    /// marks: the clone of each cell marked at the tracing step, counted as
    /// an experiment counts basal cells, cycling, quiescent, or committed to
    /// differentiate and not yet gone. Empty before marking and without
    /// tracing.
    pub fn traced_clones(&self) -> Vec<(u32, usize)> {
        let mut counts = std::collections::BTreeMap::new();
        for (i, mark) in self.marks.iter().enumerate() {
            let Some(mark) = mark else { continue };
            if self.in_layer(cell_of(i)) {
                *counts.entry(*mark).or_insert(0) += 1;
            }
        }
        counts.into_iter().collect()
    }

    /// The mean level of the evasion node over the living tumour cells, or
    /// `None` when their networks name none or none is living.
    pub fn evasion(&self) -> Option<f64> {
        let levels: Vec<f64> = self
            .networks
            .iter()
            .enumerate()
            .filter(|(i, _)| {
                let id = cell_of(*i);
                self.lattice.cell_kind(id) == Some(CellKind::Tumor)
                    && matches!(self.fate(id), Some(Fate::Cycling | Fate::Quiescent))
            })
            .filter_map(|(_, slot)| slot.as_ref().and_then(Network::evasion))
            .collect();
        (!levels.is_empty()).then(|| levels.iter().sum::<f64>() / levels.len() as f64)
    }

    /// The chemokine living tumour cells secrete, as it was last settled,
    /// when immune cells come for them.
    pub fn chemokine(&self) -> Option<&ScalarField> {
        self.chemokine.as_ref()
    }

    /// Cells a growing vessel left with no pixel at all, which died of it.
    /// Counted among the gone once resorbed, which is at once.
    pub fn crushed(&self) -> u64 {
        self.crushed
    }
}

/// What the vessel at `(x, y)` delivers, between zero and one, when what it
/// delivers follows the blood through it: nothing without blood, half at a
/// flow of `half`, as a share of what the capillary alone carries, and all
/// of it at the capillary's own flow,
/// $\frac{q (1 + h)}{q + h}$ for a flow $q$ and `half` $h$.
fn delivered(vessels: &Vessels, x: u32, y: u32, half: f64) -> f64 {
    if !vessels.carries_blood(x, y) {
        return 0.0;
    }
    let q = vessels.flow_through(x, y);
    (q * (1.0 + half) / (q + half)).min(1.0)
}

/// Refuse what a volume cannot run yet: what is built on vessels laid in the
/// plane, on a sheet, or on places named in the plane. A tumour in a bath or
/// around a straight vessel, mutating or exposed to a stress, runs in either.
fn refuse_in_a_volume(scenario: &Scenario) -> Result<(), CellweaveError> {
    let refused = if matches!(scenario.oxygen_source, OxygenSource::Capillary { .. }) {
        Some("a capillary lying in the plane")
    } else if matches!(scenario.seeding, Seeding::Sheet { .. }) {
        Some("a sheet, which is a layer one cell thick")
    } else if scenario.angiogenesis.is_some() {
        Some("vessels that grow")
    } else if scenario.immunity.is_some() {
        Some("immune cells")
    } else if scenario.lesion.is_some() {
        Some("a lesion")
    } else {
        None
    };
    match refused {
        Some(what) => Err(CellweaveError::Config(format!(
            "{what} runs only in the plane for now; set depth = 1"
        ))),
        None => Ok(()),
    }
}

/// The row a capillary lying across the box runs along: a quarter of the way
/// down, so that the tissue growing below it has most of the box.
pub fn capillary_row(height: u32) -> u32 {
    height / 4
}

/// The decay rate that gives VEGF the reach a scenario states.
fn vegf_decay(angiogenesis: &Angiogenesis) -> f64 {
    VEGF_DIFFUSION / (angiogenesis.reach * angiogenesis.reach)
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

/// What a cell asks to be at step `mcs` of a cycle that began at `began`
/// and ends at `due`: half of `full` at birth, rising in proportion to the
/// time passed, all of it from `due` on. A daughter is born at about half
/// her mother's volume, and grows through her cycle to divide again. A
/// cycle that is behind the cell, `due` not after `began`, gives all of it.
fn grown(full: Volume, began: u64, due: u64, mcs: u64) -> Volume {
    if due <= began {
        return full;
    }
    let passed = (mcs.saturating_sub(began) as f64 / (due - began) as f64).min(1.0);
    Volume((f64::from(full.0) * (0.5 + 0.5 * passed)).round() as u32)
}

/// Settle a field that every pixel of a cell of the kinds `releases` names
/// releases `production` into, a step's worth, and that the medium the bath
/// reaches holds at zero: the waste, from the living, and the necrotic
/// factor, from the dead. `cells` is how many cells the tissue has had.
fn settle_released(
    field: &mut ScalarField,
    lattice: &CpmLattice,
    cells: usize,
    production: f64,
    releases: impl Fn(CellId, CellKind) -> bool,
) -> Result<(), CellweaveError> {
    let mut per_cell = vec![0.0; cells + 1];
    for (id, record) in lattice.cell_records() {
        if releases(id, record.kind)
            && let Some(rate) = per_cell.get_mut(id.0 as usize)
        {
            *rate = production;
        }
    }
    apply_production(field, lattice, &per_cell);
    field.release_all();
    hold_medium_bath(field, lattice, 0.0);
    field.relax(
        TimeStep(1.0),
        OXYGEN_TOLERANCE * production,
        OXYGEN_MAX_STEPS,
    )?;
    Ok(())
}

/// The cell whose slot is `i`. Identifier 0 is the medium.
fn cell_of(i: usize) -> CellId {
    CellId(i as u32 + 1)
}

#[cfg(test)]
mod tests {
    use super::*;
    use cellweave_core::InTheWay;

    fn spec() -> NetworkSpec {
        toml::from_str(include_str!("../../../examples/networks/ras-erk.toml"))
            .expect("the shipped network parses")
    }

    fn scenario(oxygen_source: OxygenSource) -> Scenario {
        Scenario {
            width: 80,
            height: 80,
            depth: 1,
            max_cells: DEFAULT_MAX_CELLS,
            temperature: 10.0,
            oxygen_source,
            oxygen_uptake: 0.0003,
            seed: DEFAULT_SEED,
            mutation_rate: 0.0,
            suppress: Vec::new(),
            suppress_before: None,
            rewiring: 0.0,
            seeding: Seeding::Single,
            exposure: None,
            lesion: None,
            revert: None,
            angiogenesis: None,
            immunity: None,
            tracing: None,
            debris_loss: None,
            sweeps_per_step: 1,
            cycle: None,
            contact_order: 1,
            waste: None,
            necrotic_factor: None,
            still_medium: None,
        }
    }

    /// A tumour in a bath, in a box 40 voxels on every side.
    fn in_a_volume(oxygen_source: OxygenSource) -> Scenario {
        Scenario {
            width: 40,
            height: 40,
            depth: 40,
            ..scenario(oxygen_source)
        }
    }

    #[test]
    fn in_a_volume_the_first_cell_is_a_ball_at_the_centre() {
        let sim = Simulation::new(in_a_volume(OxygenSource::Bath), &spec()).unwrap();
        assert_eq!(sim.lattice().volume(CellId(1)), Some(Volume(123)));
        // Its centre voxel, in the middle slice.
        assert_eq!(sim.lattice().occupant_at((20 * 40 + 20) * 40 + 20), 1);
        let m = sim.measure();
        // A ball of 123 voxels has the radius of a sphere of that volume.
        assert!((m.tissue_radius - 3.08).abs() < 0.01, "{}", m.tissue_radius);
    }

    #[test]
    fn a_tumour_in_a_volume_grows_and_one_seed_gives_one_run() {
        let run = || {
            let (mut sim, _) = Simulation::new(in_a_volume(OxygenSource::Bath), &spec())
                .unwrap()
                .start()
                .unwrap();
            for _ in 0..120 {
                sim.step().unwrap();
            }
            (sim.living(), format!("{:?}", sim.snapshot()))
        };
        let (living, first) = run();
        assert!(living > 1, "nothing divided");
        assert_eq!(first, run().1);
    }

    #[test]
    fn a_still_medium_sits_below_the_source_by_what_the_tissue_takes_up() {
        let with = |still_medium| Scenario {
            oxygen_uptake: 0.01,
            still_medium,
            ..in_a_volume(OxygenSource::Bath)
        };
        let start = |scenario| {
            Simulation::new(scenario, &spec())
                .unwrap()
                .start()
                .unwrap()
                .0
        };
        assert_eq!(start(with(None)).bath_level(), 1.0);

        // The first cell, a ball of 123 voxels each taking up 0.01, in a
        // medium through which oxygen diffuses as through the tissue.
        let mut sim = start(with(Some(StillMedium { diffusion: 0.2 })));
        let radius = (3.0 * 123.0 / (4.0 * std::f64::consts::PI)).cbrt();
        let expected = 1.0 - 0.01 * 123.0 / (4.0 * std::f64::consts::PI * 0.2 * radius);
        let level = sim.bath_level();
        assert!(
            (level - expected).abs() < 1e-6,
            "{level} against {expected}"
        );
        // The step holds the medium at it: the box's corner is medium.
        sim.step().unwrap();
        assert_eq!(sim.oxygen().value_at(0), level);

        // A medium oxygen goes through faster is depleted less.
        let faster = start(with(Some(StillMedium { diffusion: 2.0 }))).bath_level();
        assert!(faster > level && faster < 1.0, "{faster} against {level}");
    }

    #[test]
    fn a_still_medium_is_a_bath_around_a_ball_that_lets_oxygen_through() {
        let still = Some(StillMedium { diffusion: 0.3 });
        let refused = |scenario, what: &str| match Simulation::new(scenario, &spec()) {
            Err(CellweaveError::Config(m)) => assert!(m.contains(what), "{m}"),
            Err(e) => panic!("{e}"),
            Ok(_) => panic!("accepted, where it should say {what}"),
        };
        refused(
            Scenario {
                still_medium: still,
                ..scenario(OxygenSource::Bath)
            },
            "a plane has none",
        );
        refused(
            Scenario {
                still_medium: still,
                ..in_a_volume(OxygenSource::Vessel { radius: 4 })
            },
            "the bath's",
        );
        refused(
            Scenario {
                still_medium: Some(StillMedium { diffusion: 0.0 }),
                ..in_a_volume(OxygenSource::Bath)
            },
            "above zero",
        );
        refused(
            Scenario {
                still_medium: still,
                waste: Some(Waste { production: 0.001 }),
                ..in_a_volume(OxygenSource::Bath)
            },
            "only oxygen",
        );
    }

    #[test]
    fn around_a_vessel_in_a_volume_the_vessel_runs_through_the_box() {
        let sim =
            Simulation::new(in_a_volume(OxygenSource::Vessel { radius: 4 }), &spec()).unwrap();
        let lattice = sim.lattice();
        let vessel = lattice.grid.iter().filter(|&&v| v == OBSTACLE).count();
        let per_slice = (0..1600)
            .filter(|&i| lattice.occupant_at(i) == OBSTACLE)
            .count();
        assert_eq!(vessel, 40 * per_slice);
        // The first cell sits against its wall, in the middle slice.
        assert_eq!(lattice.occupant_at((20 * 40 + 20) * 40 + 20 + 4 + 3), 1);
    }

    #[test]
    fn what_runs_only_in_the_plane_is_refused_in_a_volume() {
        let refused = |scenario: Scenario| {
            matches!(
                Simulation::new(scenario, &spec()),
                Err(CellweaveError::Config(m)) if m.contains("only in the plane")
            )
        };
        assert!(refused(in_a_volume(OxygenSource::Capillary { radius: 3 })));
        assert!(refused(Scenario {
            seeding: Seeding::Sheet { side: 10 },
            oxygen_uptake: 0.0,
            ..in_a_volume(OxygenSource::Bath)
        }));
        assert!(refused(Scenario {
            angiogenesis: growing(0.3).angiogenesis,
            ..in_a_volume(OxygenSource::Vessel { radius: 4 })
        }));
        assert!(refused(Scenario {
            immunity: immune(0, 0.5, 0.5).immunity,
            ..in_a_volume(OxygenSource::Vessel { radius: 4 })
        }));
    }

    #[test]
    fn a_dead_cell_loses_the_share_of_its_volume_debris_loss_says() {
        let loss = 0.05;
        let (mut sim, _) = Simulation::new(
            Scenario {
                debris_loss: Some(loss),
                ..scenario(OxygenSource::Bath)
            },
            &spec(),
        )
        .unwrap()
        .start()
        .unwrap();
        // Grown first, so that the cell at the centre has living cells all
        // around it to close over what it loses.
        for _ in 0..100 {
            sim.step().unwrap();
        }
        let id = CellId(sim.lattice().occupant(40, 40));
        assert!(id.0 != 0);
        sim.lattice.set_kind(id, CellKind::Necrotic);
        let before = f64::from(sim.lattice().volume(id).unwrap().0);
        for _ in 0..10 {
            sim.step().unwrap();
        }
        let after = f64::from(sim.lattice().volume(id).unwrap().0);
        let expected = before * (1.0 - loss).powi(10);
        assert!(
            (after - expected).abs() < 0.15 * expected,
            "{before} pixels became {after}, against {expected:.0}"
        );
    }

    #[test]
    fn with_a_debris_loss_a_cell_is_kept_from_the_step_it_dies() {
        // Ten times the uptake, so that the middle of the disc dies.
        let (mut sim, _) = Simulation::new(
            Scenario {
                debris_loss: Some(0.05),
                oxygen_uptake: 0.003,
                ..scenario(OxygenSource::Bath)
            },
            &spec(),
        )
        .unwrap()
        .start()
        .unwrap();
        let mut dead = 0;
        for _ in 0..300 {
            let mcs = sim.step().unwrap().mcs;
            for (i, slot) in sim.networks.iter().enumerate() {
                let id = cell_of(i);
                if slot.is_some()
                    && sim.lattice.cell_kind(id) == Some(CellKind::Necrotic)
                    && sim.lattice.volume(id).is_some_and(|v| v.0 > 0)
                {
                    dead += 1;
                    assert!(
                        sim.kept_dead.contains(&i),
                        "cell {i} dead at MCS {mcs}, not kept"
                    );
                }
            }
            if dead > 0 {
                break;
            }
        }
        assert!(dead > 0, "no cell died");
    }

    #[test]
    fn a_cell_grows_through_its_cycle_from_half_its_target_to_all_of_it() {
        let full = Volume(180);
        assert_eq!(grown(full, 100, 120, 100), Volume(90));
        assert_eq!(grown(full, 100, 120, 110), Volume(135));
        assert_eq!(grown(full, 100, 120, 120), full);
        assert_eq!(grown(full, 100, 120, 150), full, "arrested past its cycle");
        // A cell the tissue starts with has its cycle behind it.
        assert_eq!(grown(full, 0, 0, 7), full);
    }

    #[test]
    fn a_cycle_lasts_a_step_at_least_and_spreads_by_something() {
        for (median, spread) in [
            (0.5, 0.2),
            (20.0, -0.1),
            (f64::NAN, 0.2),
            (20.0, f64::INFINITY),
        ] {
            let scenario = Scenario {
                cycle: Some(Cycle { median, spread }),
                ..scenario(OxygenSource::Bath)
            };
            assert!(
                Simulation::new(scenario, &spec()).is_err(),
                "{median} {spread}"
            );
        }
    }

    #[test]
    fn cycle_lengths_have_the_median_and_spread_asked_for() {
        let cycle = Cycle {
            median: 40.0,
            spread: 0.3,
        };
        let (mut sim, _) = Simulation::new(scenario(OxygenSource::Bath), &spec())
            .unwrap()
            .start()
            .unwrap();
        let mut logs: Vec<f64> = (0..4000)
            .map(|_| (sim.cycle_length(cycle) as f64).ln())
            .collect();
        logs.sort_by(f64::total_cmp);
        let median = logs[logs.len() / 2].exp();
        let mean = logs.iter().sum::<f64>() / logs.len() as f64;
        let sd = (logs.iter().map(|l| (l - mean).powi(2)).sum::<f64>() / logs.len() as f64).sqrt();
        assert!((median - 40.0).abs() < 2.0, "median {median}");
        assert!((sd - 0.3).abs() < 0.03, "spread {sd}");
    }

    #[test]
    fn with_a_cycle_of_their_own_cells_stop_dividing_together() {
        // Without a cycle, cells in a volume regrow within a few steps and
        // all divide at the same checks, every 20 steps; with one, each
        // divides on its own clock.
        let steps_with_births = |cycle: Option<Cycle>| {
            let (mut sim, _) = Simulation::new(
                Scenario {
                    cycle,
                    ..in_a_volume(OxygenSource::Bath)
                },
                &spec(),
            )
            .unwrap()
            .start()
            .unwrap();
            let mut at = BTreeSet::new();
            let mut before = sim.born();
            for _ in 0..160 {
                let mcs = sim.step().unwrap().mcs;
                if sim.born() > before {
                    at.insert(mcs % DIVISION_CHECK_EVERY);
                }
                before = sim.born();
            }
            at.len()
        };
        assert_eq!(steps_with_births(None), 1);
        let cycle = Cycle {
            median: 20.0,
            spread: 0.3,
        };
        assert!(steps_with_births(Some(cycle)) >= 10);
    }

    #[test]
    fn the_contact_order_is_one_two_or_three() {
        for order in [0, 4] {
            let scenario = Scenario {
                contact_order: order,
                ..scenario(OxygenSource::Bath)
            };
            assert!(Simulation::new(scenario, &spec()).is_err(), "{order}");
        }
    }

    #[test]
    fn a_waste_is_released_at_some_rate_and_washed_away_by_a_bath() {
        for production in [0.0, -1.0, f64::NAN] {
            let scenario = Scenario {
                waste: Some(Waste { production }),
                ..scenario(OxygenSource::Bath)
            };
            assert!(Simulation::new(scenario, &spec()).is_err(), "{production}");
        }
        let scenario = Scenario {
            waste: Some(Waste { production: 0.001 }),
            ..scenario(OxygenSource::Vessel { radius: 3 })
        };
        assert!(Simulation::new(scenario, &spec()).is_err());
    }

    #[test]
    fn the_waste_is_deepest_where_the_tissue_is_and_nothing_in_the_bath() {
        let (mut sim, _) = Simulation::new(
            Scenario {
                waste: Some(Waste { production: 0.001 }),
                ..scenario(OxygenSource::Bath)
            },
            &spec(),
        )
        .unwrap()
        .start()
        .unwrap();
        for _ in 0..50 {
            sim.step().unwrap();
        }
        let field = sim.waste().unwrap();
        let at = |x: usize, y: usize| field.value_at(y * 80 + x);
        assert_eq!(sim.lattice().occupant(0, 0), 0);
        assert_eq!(at(0, 0), 0.0, "the bath washes it away");
        // Out from the middle of the tissue along a row, to the medium.
        let edge = (40..80)
            .find(|&x| sim.lattice().occupant(x as u32, 40) == 0)
            .unwrap();
        assert!(edge > 44, "the tissue reaches {edge}");
        let (middle, halfway) = (at(40, 40), at((40 + edge) / 2, 40));
        assert!(
            middle > halfway && halfway > 0.0 && at(79, 40) == 0.0,
            "{middle} at the middle, {halfway} halfway, {} at the box's edge",
            at(79, 40)
        );
    }

    #[test]
    fn a_necrotic_factor_is_released_at_some_rate_and_washed_away_by_a_bath() {
        for production in [0.0, -1.0, f64::NAN] {
            let scenario = Scenario {
                necrotic_factor: Some(NecroticFactor {
                    production,
                    release_for: None,
                }),
                ..scenario(OxygenSource::Bath)
            };
            assert!(Simulation::new(scenario, &spec()).is_err(), "{production}");
        }
        let scenario = Scenario {
            necrotic_factor: Some(NecroticFactor {
                production: 0.001,
                release_for: None,
            }),
            ..scenario(OxygenSource::Vessel { radius: 3 })
        };
        assert!(Simulation::new(scenario, &spec()).is_err());
    }

    #[test]
    fn a_dead_cell_releases_the_necrotic_factor_for_as_many_steps_as_it_is_told() {
        let scenario = |release_for| Scenario {
            necrotic_factor: Some(NecroticFactor {
                production: 0.001,
                release_for,
            }),
            ..scenario(OxygenSource::Bath)
        };
        assert!(Simulation::new(scenario(Some(0)), &spec()).is_err());
        // The same tissue either way, since nothing reads the factor; the
        // centre's level at each of six steps from the one a cell is seen dead.
        let centre = |release_for| {
            let (mut sim, _) = Simulation::new(scenario(release_for), &spec())
                .unwrap()
                .start()
                .unwrap();
            for _ in 0..50 {
                sim.step().unwrap();
            }
            let id = CellId(sim.lattice().occupant(40, 40));
            sim.lattice.set_kind(id, CellKind::Necrotic);
            (0..6)
                .map(|_| {
                    sim.step().unwrap();
                    sim.necrotic_factor().unwrap().value_at(40 * 80 + 40)
                })
                .collect::<Vec<_>>()
        };
        let (limited, always) = (centre(Some(2)), centre(None));
        assert_eq!(
            limited[..3],
            always[..3],
            "released at the step it is seen dead and two more"
        );
        assert!(
            limited[5] < 0.5 * always[5],
            "then no longer: {limited:?} against {always:?}"
        );
    }

    #[test]
    fn the_necrotic_factor_comes_from_the_dead_alone() {
        let (mut sim, _) = Simulation::new(
            Scenario {
                necrotic_factor: Some(NecroticFactor {
                    production: 0.001,
                    release_for: None,
                }),
                ..scenario(OxygenSource::Bath)
            },
            &spec(),
        )
        .unwrap()
        .start()
        .unwrap();
        for _ in 0..50 {
            sim.step().unwrap();
        }
        let field = sim.necrotic_factor().unwrap();
        assert!(
            (0..field.len()).all(|i| field.value_at(i) == 0.0),
            "nothing is dead yet, and nothing is released"
        );
        let id = CellId(sim.lattice().occupant(40, 40));
        sim.lattice.set_kind(id, CellKind::Necrotic);
        sim.step().unwrap();
        let field = sim.necrotic_factor().unwrap();
        let at = |x: usize, y: usize| field.value_at(y * 80 + x);
        let edge = (40..80)
            .find(|&x| sim.lattice().occupant(x as u32, 40) == 0)
            .unwrap();
        let (dead, halfway) = (at(40, 40), at((40 + edge) / 2, 40));
        assert!(
            dead > halfway && halfway > 0.0 && at(0, 0) == 0.0,
            "{dead} in the dead cell, {halfway} halfway out, {} in the bath",
            at(0, 0)
        );
    }

    #[test]
    fn a_step_sweeps_the_lattice_at_least_once() {
        let scenario = Scenario {
            sweeps_per_step: 0,
            ..scenario(OxygenSource::Bath)
        };
        assert!(Simulation::new(scenario, &spec()).is_err());
    }

    #[test]
    fn a_debris_loss_is_a_share_per_step() {
        for loss in [0.0, -0.1, 1.5, f64::NAN] {
            let scenario = Scenario {
                debris_loss: Some(loss),
                ..scenario(OxygenSource::Bath)
            };
            assert!(Simulation::new(scenario, &spec()).is_err(), "{loss}");
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

    #[test]
    fn undoing_a_neutral_mutation_leaves_the_history_identical_to_the_last_bit() {
        let (real, grid, clones) = history(mutating(), 160);
        let neutral = clones
            .iter()
            .find(|c| c.founded_at < 120 && c.mutation.as_ref().is_some_and(|m| m.node == "vegf"))
            .expect("one of the early clones changed VEGF");
        let revert = Some(Reversion {
            clones: vec![neutral.id],
            at: 120,
        });
        let (undone, other_grid, records) = history(
            Scenario {
                revert,
                ..mutating()
            },
            160,
        );
        assert_eq!(real, undone);
        assert_eq!(grid, other_grid);
        assert_eq!(records[neutral.id.0 as usize].reverted_at, Some(120));
    }

    #[test]
    fn undoing_a_mutation_parts_the_histories_only_from_its_step() {
        let (real, _, clones) = history(mutating(), 160);
        // Every clone founded before step 120 undone at once: something they
        // did has to show.
        let early: Vec<CloneId> = clones
            .iter()
            .filter(|c| c.mutation.is_some() && c.founded_at < 120)
            .map(|c| c.id)
            .collect();
        assert!(early.len() > 2);
        let (undone, _, _) = history(
            Scenario {
                revert: Some(Reversion {
                    clones: early,
                    at: 120,
                }),
                ..mutating()
            },
            160,
        );
        assert_eq!(
            real[..120],
            undone[..120],
            "the same history until step 120"
        );
        assert_ne!(real, undone, "and a different one after");
    }

    #[test]
    fn every_mutation_before_a_step_is_drawn_and_not_made() {
        let (real, _, clones) = history(mutating(), 100);
        let before = 60;
        let (without, _, kept) = history(
            Scenario {
                suppress_before: Some(before),
                ..mutating()
            },
            100,
        );
        assert!(
            kept.iter().any(|c| c.founded_at < before)
                && kept.iter().any(|c| c.founded_at >= before)
        );
        for clone in kept.iter().skip(1) {
            assert_eq!(clone.suppressed, clone.founded_at < before, "{clone:?}");
        }
        let first = clones[1].founded_at as usize;
        assert_eq!(
            real[..=first],
            without[..=first],
            "one history until the first mutation"
        );
    }

    #[test]
    fn coverage_counts_the_living_whose_lineage_carries_a_change() {
        let (mut sim, _) = Simulation::new(mutating(), &spec())
            .unwrap()
            .start()
            .unwrap();
        for _ in 0..100 {
            sim.step().unwrap();
        }
        let sizes = sim.clone_sizes();
        let living: usize = sizes.iter().sum();
        let coverage = sim.coverage();
        assert_eq!(coverage.len(), 5, "one share per node");
        assert!(coverage.iter().any(|(_, share)| *share > 0.0));
        for (node, share) in coverage {
            let carries = |c: usize| {
                let mut at = sim.clones().get(c);
                while let Some(record) = at {
                    if record.mutation.as_ref().is_some_and(|m| m.node == node) {
                        return true;
                    }
                    at = record.parent.and_then(|p| sim.clones().get(p.0 as usize));
                }
                false
            };
            let carrying: usize = (0..sizes.len())
                .filter(|&c| carries(c))
                .map(|c| sizes[c])
                .sum();
            assert_eq!(share, carrying as f64 / living as f64, "{node}");
        }
    }

    #[test]
    fn a_suppressed_mutation_is_not_carried() {
        let scenario = Scenario {
            suppress: (1..1000).map(CloneId).collect(),
            ..mutating()
        };
        let (mut sim, _) = Simulation::new(scenario, &spec()).unwrap().start().unwrap();
        for _ in 0..100 {
            sim.step().unwrap();
        }
        assert!(sim.clones().len() > 3, "clones were founded");
        assert!(sim.coverage().iter().all(|(_, share)| *share == 0.0));
    }

    #[test]
    fn a_ceiling_below_the_division_ratio_is_refused() {
        let mut spec = spec();
        if let Some(m) = spec.mechanics.as_mut() {
            m.max_target = Some(1.2);
        }
        assert!(Simulation::new(scenario(OxygenSource::Bath), &spec).is_err());
    }

    #[test]
    fn with_full_rewiring_every_mutation_gains_or_loses_a_term() {
        let (_, _, clones) = history(
            Scenario {
                rewiring: 1.0,
                ..mutating()
            },
            100,
        );
        assert!(clones.len() > 3, "the population grew and mutated");
        for clone in clones.iter().skip(1) {
            let change = &clone.mutation.as_ref().unwrap().change;
            assert!(
                !matches!(change, crate::network::Change::Scaled(_)),
                "{change:?}"
            );
        }
    }

    fn growing(threshold: f64) -> Scenario {
        Scenario {
            angiogenesis: Some(Angiogenesis {
                threshold,
                speed: 1.0,
                reach: 15.0,
                perfusion: Perfusion::AtOnce,
                regress_after: None,
                half_delivery_flow: None,
                prune_below: None,
                in_the_way: InTheWay::Crushed,
            }),
            width: 100,
            height: 100,
            oxygen_uptake: 0.0001,
            ..scenario(OxygenSource::Vessel { radius: 4 })
        }
    }

    fn vessel_pixels(sim: &Simulation<Running>) -> usize {
        sim.lattice()
            .grid
            .iter()
            .filter(|&&p| p == OBSTACLE)
            .count()
    }

    #[test]
    fn a_vessel_that_never_sprouts_leaves_the_run_as_it_was() {
        let run = |scenario: Scenario| {
            let (mut sim, _) = Simulation::new(scenario, &spec()).unwrap().start().unwrap();
            for _ in 0..80 {
                sim.step().unwrap();
            }
            (format!("{:?}", sim.snapshot()), sim.lattice().grid.clone())
        };
        assert_eq!(
            run(growing(1e9)),
            run(Scenario {
                angiogenesis: None,
                immunity: None,
                tracing: None,
                ..growing(1e9)
            })
        );
    }

    #[test]
    fn starving_tissue_grows_vessel_towards_itself() {
        let (mut sim, _) = Simulation::new(growing(0.3), &spec())
            .unwrap()
            .start()
            .unwrap();
        let before = vessel_pixels(&sim);
        for _ in 0..300 {
            if sim.step().unwrap().reached_edge {
                break;
            }
        }
        assert!(sim.sprouted() > 0, "something called and the wall sprouted");
        assert!(vessel_pixels(&sim) > before + 30, "and laid vessel");
        assert!(sim.vegf().is_some());
    }

    #[test]
    fn vessels_grow_only_from_a_vessel_towards_a_network_that_calls() {
        let refused = |scenario: Scenario, spec: &NetworkSpec| {
            matches!(
                Simulation::new(scenario, spec),
                Err(CellweaveError::Config(_))
            )
        };
        let angiogenesis = growing(0.3).angiogenesis;
        assert!(refused(
            Scenario {
                angiogenesis,
                ..scenario(OxygenSource::Bath)
            },
            &spec()
        ));
        let mut silent = spec();
        if let Some(roles) = silent.roles.as_mut() {
            roles.angiogenic = None;
        }
        assert!(refused(growing(0.3), &silent));
        for (speed, reach) in [(1.5, 15.0), (0.0, 15.0), (1.0, 0.0)] {
            let mut s = growing(0.3);
            s.angiogenesis = Some(Angiogenesis {
                threshold: 0.3,
                speed,
                reach,
                perfusion: Perfusion::AtOnce,
                regress_after: None,
                half_delivery_flow: None,
                prune_below: None,
                in_the_way: InTheWay::Crushed,
            });
            assert!(refused(s, &spec()), "speed {speed}, reach {reach}");
        }
        assert!(refused(growing(0.0), &spec()));
    }

    fn dry_for(regress_after: Option<u64>) -> Scenario {
        let mut s = growing(0.3);
        s.angiogenesis = s.angiogenesis.map(|a| Angiogenesis {
            perfusion: Perfusion::OnJoining,
            regress_after,
            ..a
        });
        s
    }

    #[test]
    fn a_sprout_regresses_only_where_it_can_go_without_blood_for_a_while() {
        let refused = |scenario: Scenario| {
            matches!(
                Simulation::new(scenario, &spec()),
                Err(CellweaveError::Config(_))
            )
        };
        assert!(refused(dry_for(Some(0))));
        let mut at_once = dry_for(Some(30));
        at_once.angiogenesis = at_once.angiogenesis.map(|a| Angiogenesis {
            perfusion: Perfusion::AtOnce,
            ..a
        });
        assert!(refused(at_once));
        assert!(!refused(dry_for(Some(30))));
    }

    #[test]
    fn tissue_fills_where_a_blind_sprout_regressed() {
        let (mut sim, _) = Simulation::new(dry_for(Some(10)), &spec())
            .unwrap()
            .start()
            .unwrap();
        let laid_at_start: Vec<bool> = sim.lattice().grid.iter().map(|&p| p == OBSTACLE).collect();
        let mut ever_sprouted = vec![false; laid_at_start.len()];
        for _ in 0..300 {
            if sim.step().unwrap().reached_edge {
                break;
            }
            for (seen, (&now, &first)) in ever_sprouted
                .iter_mut()
                .zip(sim.lattice().grid.iter().zip(&laid_at_start))
            {
                *seen |= now == OBSTACLE && !first;
            }
        }
        assert!(sim.regressed() > 0, "a blind sprout regressed");
        let refilled = sim
            .lattice()
            .grid
            .iter()
            .zip(&ever_sprouted)
            .filter(|&(&now, &seen)| seen && now != 0 && now != OBSTACLE)
            .count();
        assert!(refilled > 0, "and cells moved in where it was");
    }

    #[test]
    fn a_capillary_lies_across_the_box_and_the_tissue_grows_below_it() {
        let (mut sim, _) =
            Simulation::new(scenario(OxygenSource::Capillary { radius: 2 }), &spec())
                .unwrap()
                .start()
                .unwrap();
        let row = capillary_row(80);
        let lattice = sim.lattice();
        assert_eq!(lattice.occupant(0, row), OBSTACLE, "from the left side");
        assert_eq!(lattice.occupant(79, row), OBSTACLE, "to the right side");
        assert_ne!(
            lattice.occupant(40, row + 2 + 8),
            0,
            "the first cell below it"
        );
        for _ in 0..100 {
            sim.step().unwrap();
        }
        let above = (0..row.saturating_sub(2))
            .flat_map(|y| (0..80).map(move |x| (x, y)))
            .any(|(x, y)| {
                let id = sim.lattice().occupant(x, y);
                id != 0 && id != OBSTACLE
            });
        assert!(!above, "nothing crosses the capillary");
        let m = sim.measure();
        assert!(m.vessel > 0 && m.viable_depth > 0.0, "{m:?}");
    }

    #[test]
    fn blood_flows_as_pressure_drives_it_only_through_a_capillary() {
        let flowing = |oxygen_source: OxygenSource| {
            let mut s = growing(0.3);
            s.oxygen_source = oxygen_source;
            s.angiogenesis = s.angiogenesis.map(|a| Angiogenesis {
                perfusion: Perfusion::Flow,
                ..a
            });
            s
        };
        assert!(matches!(
            Simulation::new(flowing(OxygenSource::Vessel { radius: 4 }), &spec()),
            Err(CellweaveError::Config(_))
        ));
        let (mut sim, _) = Simulation::new(flowing(OxygenSource::Capillary { radius: 2 }), &spec())
            .unwrap()
            .start()
            .unwrap();
        for _ in 0..200 {
            if sim.step().unwrap().reached_edge {
                break;
            }
        }
        let m = sim.measure();
        assert!(m.perfused <= m.vessel);
        let row = capillary_row(100);
        assert!(sim.carries_blood(50, row), "the capillary always does");
    }

    #[test]
    fn delivery_by_flow_needs_flow_and_a_half_above_zero() {
        let delivering = |perfusion: Perfusion, half: f64| {
            let mut s = growing(0.3);
            s.oxygen_source = OxygenSource::Capillary { radius: 2 };
            s.angiogenesis = s.angiogenesis.map(|a| Angiogenesis {
                perfusion,
                half_delivery_flow: Some(half),
                ..a
            });
            Simulation::new(s, &spec())
        };
        assert!(delivering(Perfusion::OnJoining, 0.05).is_err());
        assert!(delivering(Perfusion::Flow, 0.0).is_err());
        // Before anything grows, the capillary carries its own flow and
        // delivers in full.
        let (sim, _) = delivering(Perfusion::Flow, 0.05).unwrap().start().unwrap();
        let row = capillary_row(100);
        let at = |y| cellweave_core::DiffusionField::concentration_at(sim.oxygen(), 50, y);
        assert_eq!(at(row), 1.0);
        assert!(at(row + 12) > 0.5, "and feeds the tissue below");
    }

    #[test]
    fn pruning_by_flow_needs_flow_regression_and_a_floor_above_zero() {
        let pruning = |perfusion: Perfusion, regress_after: Option<u64>, floor: f64| {
            let mut s = growing(0.3);
            s.oxygen_source = OxygenSource::Capillary { radius: 2 };
            s.angiogenesis = s.angiogenesis.map(|a| Angiogenesis {
                perfusion,
                regress_after,
                prune_below: Some(floor),
                ..a
            });
            Simulation::new(s, &spec())
        };
        assert!(pruning(Perfusion::OnJoining, Some(100), 0.02).is_err());
        assert!(pruning(Perfusion::Flow, None, 0.02).is_err());
        assert!(pruning(Perfusion::Flow, Some(100), 0.0).is_err());
        assert!(pruning(Perfusion::Flow, Some(100), 0.02).is_ok());
    }

    fn immune(from: u64, arrival: f64, kill: f64) -> Scenario {
        Scenario {
            immunity: Some(Immunity {
                from,
                arrival,
                kill,
                kills: 3,
                lifetime: 100,
                chemotaxis: 200.0,
                reach: 10.0,
            }),
            ..scenario(OxygenSource::Vessel { radius: 4 })
        }
    }

    #[test]
    fn immune_cells_come_only_through_a_vessel_and_only_from_their_step() {
        let refused =
            |s: Scenario| matches!(Simulation::new(s, &spec()), Err(CellweaveError::Config(_)));
        let mut bath = immune(0, 0.5, 0.5);
        bath.oxygen_source = OxygenSource::Bath;
        assert!(refused(bath));
        assert!(refused(immune(0, 0.0, 0.5)));
        assert!(refused(immune(0, 0.5, 1.5)));
        let (mut sim, _) = Simulation::new(immune(30, 1.0, 1.0), &spec())
            .unwrap()
            .start()
            .unwrap();
        for _ in 0..30 {
            sim.step().unwrap();
        }
        assert_eq!(sim.immune_cells(), 0, "none before step 30");
        sim.step().unwrap();
        assert_eq!(
            sim.immune_cells(),
            1,
            "one arrives at step 30, the chance being one"
        );
    }

    #[test]
    fn an_immune_cell_kills_the_tumour_cell_it_touches() {
        let (mut sim, _) = Simulation::new(immune(0, 1.0, 1.0), &spec())
            .unwrap()
            .start()
            .unwrap();
        for _ in 0..60 {
            sim.step().unwrap();
            if sim.killed() > 0 {
                break;
            }
        }
        assert!(
            sim.killed() > 0,
            "the first tumour cell sits by the vessel they arrive through"
        );
        let tumour = sim
            .lattice()
            .cell_records()
            .filter(|(_, r)| matches!(r.kind, CellKind::Tumor | CellKind::Normal) && r.volume.0 > 0)
            .count();
        assert_eq!(
            sim.living(),
            tumour,
            "immune cells are not counted among the living"
        );
    }

    #[test]
    fn a_tumour_cell_fully_shielded_is_never_killed() {
        let child: NetworkSpec =
            toml::from_str(include_str!("../../../examples/networks/ras-erk-pdl1.toml")).unwrap();
        let mut shielded = child.merge_onto(spec());
        if let Some(node) = shielded.nodes.iter_mut().find(|n| n.name == "pdl1") {
            node.initial = 1.0;
            node.terms.clear();
        }
        let (mut sim, _) = Simulation::new(immune(0, 1.0, 1.0), &shielded)
            .unwrap()
            .start()
            .unwrap();
        for _ in 0..60 {
            sim.step().unwrap();
        }
        assert!(sim.immune_cells() > 0, "immune cells came");
        assert_eq!(sim.killed(), 0, "and killed nothing");
        assert_eq!(sim.evasion(), Some(1.0));
    }
}
