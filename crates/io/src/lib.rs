//! cellweave-io: configuration loading and simulation output.

#![deny(unsafe_code)]

pub mod config;
pub mod output;

pub use config::SimConfig;
pub use output::JsonOutput;
