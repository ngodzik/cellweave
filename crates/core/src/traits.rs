//! Port traits: the interfaces between simulation layers.
//!
//! These traits define what each engine must do. Implementations live in
//! `cellweave-engine`. This crate has no knowledge of how they work.

use crate::error::CellweaveError;
use crate::types::{CellInputs, CellMechanics, ProteinState, SimSnapshot, TimeStep};

// ── Signaling network (ODE layer) ─────────────────────────────────────────────

/// A per-cell protein signaling network solved by an ODE integrator.
///
/// Each cell owns one instance. The network receives extracellular inputs and
/// advances its internal state by `dt` seconds.
///
/// # Example: minimal network that does nothing
/// ```
/// use cellweave_core::{SignalingNetwork, ProteinState, CellInputs, CellMechanics, TimeStep};
///
/// struct NullNetwork;
/// impl SignalingNetwork for NullNetwork {
///     fn step(&mut self, _dt: TimeStep, _inputs: &CellInputs) {}
///     fn state(&self) -> &ProteinState { todo!() }
///     fn mechanics(&self) -> CellMechanics { todo!() }
/// }
/// ```
pub trait SignalingNetwork: Send + Sync {
    /// Advance the network by `dt` given local extracellular `inputs`.
    fn step(&mut self, dt: TimeStep, inputs: &CellInputs);

    /// Current protein concentrations.
    fn state(&self) -> &ProteinState;

    /// Mechanical parameters derived from current protein state.
    ///
    /// Called by the CPM engine before each Monte Carlo step to get
    /// up-to-date adhesion and volume parameters.
    fn mechanics(&self) -> CellMechanics;
}

// ── Diffusion field (PDE layer) ───────────────────────────────────────────────

/// A scalar extracellular field solved by a finite-difference PDE integrator.
///
/// One instance per molecule (O2, VEGF, EGF, drug, …).
pub trait DiffusionField: Send + Sync {
    /// Advance the field by one time step.
    ///
    /// Cells' secretion and uptake must be applied inside this call, using
    /// the cell layout provided by the CPM engine.
    fn step(&mut self, dt: TimeStep) -> Result<(), CellweaveError>;

    /// Concentration at a given grid position.
    fn concentration_at(&self, x: u32, y: u32) -> f64;

    /// Dimensions (width, height).
    fn dims(&self) -> (u32, u32);
}

// ── Simulation output ─────────────────────────────────────────────────────────

/// Receiver of simulation snapshots. Implement to write VTK, PNG, JSON, etc.
pub trait SimOutput: Send {
    /// Called after each saved Monte Carlo step.
    fn write_snapshot(&mut self, snapshot: &SimSnapshot) -> Result<(), CellweaveError>;

    /// Called once when the simulation ends.
    fn finalize(&mut self) -> Result<(), CellweaveError> {
        Ok(())
    }
}
