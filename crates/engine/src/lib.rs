#![deny(unsafe_code)]
#![deny(missing_docs)]

//! cellweave-engine: CPM lattice, ODE integrator, PDE diffusion.

pub mod cpm;
pub mod ode;
pub mod pde;

pub use cpm::CpmLattice;
pub use ode::RasErkNetwork;
pub use pde::ScalarField;
