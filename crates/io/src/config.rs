//! TOML simulation configuration.

use cellweave_core::CellweaveError;
use serde::{Deserialize, Serialize};

/// Where the oxygen comes from, which is also what shape the tissue takes.
#[derive(Debug, Clone, Copy, PartialEq, Default, Serialize, Deserialize)]
#[serde(rename_all = "snake_case", tag = "kind")]
pub enum OxygenSource {
    /// A bath at the edge of the box, held against the surface of the tissue.
    /// One cell is seeded at the centre and grows into a spheroid, which
    /// starves from the inside once large enough. The run stops when the
    /// tissue touches the edge of the box.
    #[default]
    Bath,
    /// A vessel at the centre, held at full oxygen, that cells cannot enter.
    /// One cell is seeded against its wall and the tissue grows around it out
    /// to the depth diffusion can feed: a tumour cord. With dead cells resorbed,
    /// the cord keeps a steady thickness and turns over indefinitely.
    Vessel {
        /// Radius of the vessel, in pixels.
        radius: u32,
    },
}

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
    /// Where the oxygen comes from. Defaults to a bath when absent.
    #[serde(default)]
    pub oxygen_source: OxygenSource,
    /// What one consuming pixel takes up per unit time, as a fraction of the
    /// source level. Sets how deep tissue can be before it starves, and the
    /// right value depends on the geometry: a bathed disc is fed from its whole
    /// perimeter, a vessel from a small circle whose flux thins out with
    /// distance, so a cord needs roughly a tenth of what a spheroid does to reach
    /// the same depth.
    #[serde(default = "default_oxygen_uptake")]
    pub oxygen_uptake: f64,
}

fn default_oxygen_uptake() -> f64 {
    0.0003
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
            oxygen_source: OxygenSource::Bath,
            oxygen_uptake: default_oxygen_uptake(),
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
        assert_eq!(cfg2.oxygen_source, OxygenSource::Bath);
    }

    #[test]
    fn a_config_without_a_source_gets_the_bath() {
        // Every config written before the vessel existed keeps meaning what
        // it meant.
        let cfg = SimConfig::from_toml(
            "width = 50\nheight = 50\nmcs = 10\ntemperature = 10.0\nsnapshot_interval = 5\noutput_dir = \"o\"\n",
        )
        .unwrap();
        assert_eq!(cfg.oxygen_source, OxygenSource::Bath);
    }

    #[test]
    fn a_vessel_is_a_table_with_a_radius() {
        let cfg = SimConfig::from_toml(
            "width = 50\nheight = 50\nmcs = 10\ntemperature = 10.0\nsnapshot_interval = 5\noutput_dir = \"o\"\n\n[oxygen_source]\nkind = \"vessel\"\nradius = 6\n",
        )
        .unwrap();
        assert_eq!(cfg.oxygen_source, OxygenSource::Vessel { radius: 6 });
        assert_eq!(
            cfg.oxygen_uptake, 0.0003,
            "uptake keeps its default when absent"
        );
    }
}
