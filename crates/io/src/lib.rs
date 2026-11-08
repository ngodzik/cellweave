//! cellweave-io: configuration loading and simulation output.

#![deny(unsafe_code)]

pub mod config;
pub mod output;

pub use cellweave_core::{Exposure, OxygenSource, Schedule, Seeding};
pub use config::{LesionConfig, SimConfig, shipped_network};
pub use output::{JsonOutput, MeasurementLog, TraceLog, write_clones};
