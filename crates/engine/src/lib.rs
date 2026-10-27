#![deny(unsafe_code)]
#![deny(missing_docs)]

//! cellweave-engine: CPM lattice, a protein network read from a file and
//! integrated per cell, PDE diffusion, the coupling that carries information
//! between the lattice and the field, and the loop that steps them together.

pub mod coupling;
pub mod cpm;
pub mod network;
pub mod pde;
pub mod simulation;

pub use cpm::CpmLattice;
pub use network::Network;
pub use pde::ScalarField;
pub use simulation::{Fate, Scenario, Simulation};
