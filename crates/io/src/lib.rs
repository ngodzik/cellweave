//! cellweave-io: configuration loading and simulation output.

#![deny(unsafe_code)]

pub mod config;
pub mod output;

pub use cellweave_core::OxygenSource;
pub use config::SimConfig;
pub use output::{JsonOutput, MeasurementLog, write_clones};
