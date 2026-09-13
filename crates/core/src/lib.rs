#![deny(unsafe_code)]
#![deny(missing_docs)]

//! cellweave-core: domain types and traits.
//!
//! No I/O, no computation. Only the language of the simulation.

pub mod error;
pub mod traits;
pub mod types;

pub use error::CellweaveError;
pub use traits::*;
pub use types::*;
