//! Domain newtypes and value types.

use serde::{Deserialize, Serialize};

// ── Identity ──────────────────────────────────────────────────────────────────

/// Unique identifier for a cell.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct CellId(pub u32);

/// Identifier of a clone: every cell descended from one founder with no
/// mutation in between. The first cell founds clone 0, and each mutation
/// founds the next.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord, Serialize, Deserialize)]
pub struct CloneId(pub u32);

// ── Physical quantities ───────────────────────────────────────────────────────

/// Protein concentration (arbitrary units, normalized 0-1 unless stated).
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct Concentration(pub f64);

/// CPM Hamiltonian energy (dimensionless).
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct Energy(pub f64);

/// Simulation time step (Monte Carlo steps for CPM, seconds for ODE/PDE).
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, Serialize, Deserialize)]
pub struct TimeStep(pub f64);

/// Volume in lattice pixels.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct Volume(pub u32);

// ── Grid coordinates ──────────────────────────────────────────────────────────

/// 2D lattice position.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Pos2 {
    /// Horizontal coordinate.
    pub x: u32,
    /// Vertical coordinate.
    pub y: u32,
}

// ── Cell types ────────────────────────────────────────────────────────────────

/// Biological type of a cell, determines CPM energy parameters.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum CellKind {
    /// No cell: the background medium.
    Medium,
    /// Actively proliferating tumor cell.
    Tumor,
    /// Surrounding healthy tissue.
    Normal,
    /// Dead cell due to hypoxia or apoptosis.
    Necrotic,
    /// A cell that has committed to differentiate and is leaving the layer it
    /// divided in: in an epithelium, moving up out of the basal layer.
    Differentiated,
    /// A cell that has started its own programmed death, apoptosis, and is
    /// breaking up into fragments its neighbours clear.
    Apoptotic,
}

// ── Geometry ──────────────────────────────────────────────────────────────────

/// A lethal stress the tissue is exposed to: an abstract agent that comes in
/// from wherever oxygen does, spreads between the cells, is taken up by them,
/// and raises in each cell whatever its network says it raises.
///
/// It stands for nothing in particular. What the model can say about it is
/// how a population that varies, in a space where the stress does not reach
/// everyone equally, responds to being exposed to it.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Exposure {
    /// Level held at the source while exposure is on, as a fraction of the
    /// most the network files are written for.
    pub intensity: f64,
    /// What one living pixel takes up per unit time, as a fraction of the
    /// level at the source. The higher, the shorter the distance the stress
    /// reaches into the tissue.
    pub uptake: f64,
    /// First lattice step of exposure.
    pub start: u64,
    /// When exposure is on, from `start` on.
    pub schedule: Schedule,
}

/// When an exposure is on.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", tag = "kind")]
pub enum Schedule {
    /// On at every step from the start.
    Continuous,
    /// On for `on` steps, off for `off`, over and over.
    Periodic {
        /// Steps on, in each cycle.
        on: u64,
        /// Steps off, in each cycle.
        off: u64,
    },
}

impl Exposure {
    /// The level held at the source at lattice step `mcs`: the intensity
    /// while exposure is on, zero otherwise.
    pub fn level_at(&self, mcs: u64) -> f64 {
        if mcs < self.start {
            return 0.0;
        }
        let on = match self.schedule {
            Schedule::Continuous => true,
            Schedule::Periodic { on, off } => (mcs - self.start) % (on + off).max(1) < on,
        };
        if on { self.intensity } else { 0.0 }
    }
}

/// What the box holds at the start.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", tag = "kind")]
pub enum Seeding {
    /// One cell, placed where the oxygen source says: a tumour grown from a
    /// single cell.
    #[default]
    Single,
    /// The whole box tiled with square cells of normal tissue, each the
    /// founder of its own clone: an epithelium whose cells are all labelled at
    /// the start, as in a lineage tracing experiment.
    Sheet {
        /// Side of each square cell, in pixels.
        side: u32,
    },
}

/// Where the oxygen comes from, which is also what shape the tissue takes.
#[derive(Debug, Clone, Copy, PartialEq, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", tag = "kind")]
pub enum OxygenSource {
    /// A bath at the edge of the box, held against the surface of the tissue.
    /// One cell is seeded at the centre and grows into a spheroid, which
    /// starves from the inside once large enough. The run stops when the
    /// tissue touches the edge of the box.
    #[default]
    Bath,
    /// A vessel at the centre, held at full oxygen, that cells cannot enter.
    /// One cell is seeded against its wall and the tissue grows around it out
    /// to the depth diffusion can feed: a tumour cord. With dead cells resorbed,
    /// the cord keeps a steady thickness and turns over indefinitely.
    Vessel {
        /// Radius of the vessel, in pixels.
        radius: u32,
    },
}

// ── Protein state ─────────────────────────────────────────────────────────────

/// Snapshot of all protein concentrations inside one cell.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProteinState {
    /// Ordered concentrations matching the network's protein list.
    pub values: Vec<Concentration>,
}

impl ProteinState {
    /// Create a state initialized to zero for `n` proteins.
    pub fn zeros(n: usize) -> Self {
        Self {
            values: vec![Concentration(0.0); n],
        }
    }
}

// ── Cell inputs from environment ──────────────────────────────────────────────

/// Extracellular signals sensed by a cell at its current position.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CellInputs {
    /// Local O2 concentration.
    pub o2: Concentration,
    /// Local EGF concentration.
    pub egf: Concentration,
    /// Local VEGF concentration (secreted by neighbors).
    pub vegf: Concentration,
    /// How squeezed the cell is by its neighbours, from zero (at or above
    /// its lineage's reference size) to one (no room at all).
    pub crowding: Concentration,
    /// Local level of the lethal stress, as a fraction of what the source
    /// holds while exposure is on. Zero when there is no exposure.
    pub stress: Concentration,
}

// ── CPM parameters derived from ODE state ────────────────────────────────────

/// Mechanical parameters for one cell, derived from its protein state.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CellMechanics {
    /// Target volume in pixels.
    pub target_volume: Volume,
    /// Lambda volume (volume constraint strength).
    pub lambda_volume: f64,
    /// Cell-medium contact energy.
    pub j_medium: f64,
    /// Cell-cell contact energy (same type).
    pub j_self: f64,
}

// ── Simulation snapshot ───────────────────────────────────────────────────────

/// State of the whole simulation at one point in time, used for output.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimSnapshot {
    /// Current Monte Carlo step.
    pub mcs: u64,
    /// Per-cell data.
    pub cells: Vec<CellSnapshot>,
}

/// State of one cell at one point in time.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CellSnapshot {
    /// Cell identifier.
    pub id: CellId,
    /// Biological type at this step.
    pub kind: CellKind,
    /// Pixel count.
    pub volume: Volume,
    /// Protein concentrations.
    pub proteins: ProteinState,
    /// Mean oxygen over the cell's pixels, as a fraction of the bath.
    pub oxygen: Concentration,
    /// The clone the cell belongs to.
    pub clone: CloneId,
}
