#![deny(unsafe_code)]
#![deny(missing_docs)]

//! cellweave-engine: CPM lattice, a protein network read from a file and
//! integrated per cell, PDE diffusion, the coupling that carries information
//! between the lattice and the field, vessels that grow through the tissue,
//! and the loop that steps them together.

pub mod coupling;
pub mod cpm;
pub mod evolution;
pub mod measure;
pub mod network;
pub mod pde;
pub mod simulation;
pub mod sisters;
pub mod vessels;

pub use cpm::CpmLattice;
pub use evolution::CloneRecord;
pub use measure::{DEPTH_SHARE, Measurements};
pub use network::{Mutation, Network};
pub use pde::ScalarField;
pub use simulation::{Fate, Scenario, Simulation};
pub use sisters::SisterFates;
pub use vessels::Tip;
