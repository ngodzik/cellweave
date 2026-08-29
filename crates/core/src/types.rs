//! Domain newtypes and value types.

use serde::{Deserialize, Serialize};

// ── Identity ──────────────────────────────────────────────────────────────────

/// Unique identifier for a cell.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct CellId(pub u32);

/// Unique identifier for a protein species.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct ProteinId(pub String);

/// Unique identifier for an extracellular field (O2, VEGF, drug, …).
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct FieldId(pub String);

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

/// 3D lattice position.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Pos3 {
    /// X coordinate.
    pub x: u32,
    /// Y coordinate.
    pub y: u32,
    /// Z coordinate (depth).
    pub z: u32,
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
    /// Local drug concentrations keyed by field index.
    pub drugs: Vec<Concentration>,
}

impl CellInputs {
    /// Zero inputs, cell in isolation.
    pub fn zero(n_drugs: usize) -> Self {
        Self {
            o2: Concentration(0.0),
            egf: Concentration(0.0),
            vegf: Concentration(0.0),
            drugs: vec![Concentration(0.0); n_drugs],
        }
    }
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
}
