//! Simulation output: JSON snapshots to disk.

use cellweave_core::{CellweaveError, SimOutput, SimSnapshot};
use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;

/// Writes one JSON file per snapshot into a directory.
pub struct JsonOutput {
    dir: PathBuf,
}

impl JsonOutput {
    /// Create a new writer; `dir` is created if it does not exist.
    pub fn new(dir: impl Into<PathBuf>) -> Result<Self, CellweaveError> {
        let dir = dir.into();
        fs::create_dir_all(&dir).map_err(|e| CellweaveError::Config(format!("output dir: {e}")))?;
        Ok(Self { dir })
    }
}

impl SimOutput for JsonOutput {
    fn write_snapshot(&mut self, snapshot: &SimSnapshot) -> Result<(), CellweaveError> {
        let path = self.dir.join(format!("snapshot_{:06}.json", snapshot.mcs));
        let json = serde_json::to_string_pretty(snapshot)
            .map_err(|e| CellweaveError::Simulation(e.to_string()))?;
        let mut file =
            File::create(&path).map_err(|e| CellweaveError::Simulation(e.to_string()))?;
        file.write_all(json.as_bytes())
            .map_err(|e| CellweaveError::Simulation(e.to_string()))?;
        Ok(())
    }
}
