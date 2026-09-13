//! Domain error types.

use thiserror::Error;

/// Top-level error type for cellweave.
#[derive(Debug, Error)]
pub enum CellweaveError {
    /// Invalid configuration value.
    #[error("configuration error: {0}")]
    Config(String),

    /// Simulation reached an invalid state.
    #[error("simulation error: {0}")]
    Simulation(String),

    /// ODE solver failed to converge.
    #[error("ODE solver did not converge: {0}")]
    OdeSolver(String),

    /// Requested cell or field does not exist.
    #[error("not found: {0}")]
    NotFound(String),
}
