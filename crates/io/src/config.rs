//! TOML simulation configuration.

use cellweave_core::CellweaveError;
use serde::{Deserialize, Serialize};

/// Top-level simulation configuration (loaded from a `.toml` file).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SimConfig {
    /// Grid width in pixels.
    pub width: u32,
    /// Grid height in pixels.
    pub height: u32,
    /// Number of Monte Carlo steps to run.
    pub mcs: u64,
    /// CPM temperature.
    pub temperature: f64,
    /// How often to save a snapshot (every N MCS).
    pub snapshot_interval: u64,
    /// Output directory.
    pub output_dir: String,
}

impl SimConfig {
    /// Load config from a TOML string.
    pub fn from_toml(s: &str) -> Result<Self, CellweaveError> {
        toml::from_str(s).map_err(|e| CellweaveError::Config(e.to_string()))
    }

    /// Load config from a file path.
    pub fn from_file(path: &str) -> Result<Self, CellweaveError> {
        let s = std::fs::read_to_string(path)
            .map_err(|e| CellweaveError::Config(format!("{path}: {e}")))?;
        Self::from_toml(&s)
    }
}

impl Default for SimConfig {
    fn default() -> Self {
        Self {
            width: 200,
            height: 200,
            mcs: 1000,
            temperature: 10.0,
            snapshot_interval: 100,
            output_dir: "output".to_string(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrip_toml() {
        let cfg = SimConfig::default();
        let s = toml::to_string(&cfg).unwrap();
        let cfg2 = SimConfig::from_toml(&s).unwrap();
        assert_eq!(cfg.width, cfg2.width);
        assert_eq!(cfg.mcs, cfg2.mcs);
    }
}
