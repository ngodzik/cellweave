#![deny(unsafe_code)]
#![deny(missing_docs)]

//! cellweave-core: domain types and traits.
//!
//! No I/O, no computation. Only the language of the simulation.

pub mod types;
pub mod traits;
pub mod error;

pub use types::*;
pub use traits::*;
pub use error::CellweaveError;
