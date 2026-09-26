//! TOML simulation configuration.

use cellweave_core::{CellweaveError, NetworkSpec};
use serde::{Deserialize, Serialize};
use std::path::Path;

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
    /// File holding the protein network, relative to the configuration file.
    ///
    /// The hypothesis under test lives there, so trying another one is writing
    /// another file rather than editing and rebuilding the program. Left unset,
    /// the built-in RAS/ERK network is used, which is the same file compiled in.
    #[serde(default)]
    pub network: Option<String>,
}

/// The RAS/ERK network, compiled in so that the built-in hypothesis and one you
/// write go through exactly the same path.
pub const BUILT_IN_NETWORK: &str = include_str!("../../../examples/networks/ras-erk.toml");

fn default_oxygen_uptake() -> f64 {
    0.0003
}

impl SimConfig {
    /// Load config from a TOML string.
    pub fn from_toml(s: &str) -> Result<Self, CellweaveError> {
        toml::from_str(s).map_err(|e| CellweaveError::Config(e.to_string()))
    }

    /// Load config from a file path.
    ///
    /// # Errors
    ///
    /// [`CellweaveError::Config`] when the file cannot be read or is not a
    /// configuration.
    pub fn from_file(path: &str) -> Result<Self, CellweaveError> {
        let s = std::fs::read_to_string(path)
            .map_err(|e| CellweaveError::Config(format!("{path}: {e}")))?;
        Self::from_toml(&s)
    }

    /// The network this configuration asks for.
    ///
    /// `network` is resolved next to the configuration file, so a pair of files
    /// can be moved together. Left unset, the built-in RAS/ERK network is read
    /// from the copy compiled into the program.
    ///
    /// # Errors
    ///
    /// [`CellweaveError::Config`] when the file cannot be read or is not a
    /// network.
    pub fn load_network(&self, config_path: Option<&str>) -> Result<NetworkSpec, CellweaveError> {
        let text = match &self.network {
            None => BUILT_IN_NETWORK.to_owned(),
            Some(name) => {
                let path = match config_path.and_then(|c| Path::new(c).parent()) {
                    Some(dir) => dir.join(name),
                    None => Path::new(name).to_path_buf(),
                };
                std::fs::read_to_string(&path).map_err(|e| {
                    CellweaveError::Config(format!("reading the network {}: {e}", path.display()))
                })?
            }
        };
        toml::from_str(&text).map_err(|e| {
            let which = self.network.as_deref().unwrap_or("the built-in network");
            CellweaveError::Config(format!("{which} is not a network: {e}"))
        })
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
            network: None,
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

    #[test]
    fn the_built_in_network_is_a_network() {
        let spec = SimConfig::default().load_network(None).unwrap();
        assert_eq!(spec.nodes.len(), 5);
        assert_eq!(spec.roles.survival, "bcl2");
    }

    #[test]
    fn a_network_file_is_found_next_to_the_configuration() {
        let dir = std::env::temp_dir().join(format!("cellweave-net-{}", std::process::id()));
        std::fs::create_dir_all(dir.join("nets")).unwrap();
        std::fs::write(
            dir.join("nets").join("tiny.toml"),
            concat!(
                "[[nodes]]\nname = \"a\"\ninitial = 0.1\n\n",
                "[roles]\ngrowth = \"a\"\nsurvival = \"a\"\nhypoxia_response = \"a\"\n\n",
                "[mechanics]\ngrowth_coefficient = 1.0\nlambda_volume = 1.0\n",
                "j_medium = 1.0\nj_self_base = 1.0\nj_self_from_growth = 0.0\n",
            ),
        )
        .unwrap();
        let cfg_path = dir.join("sim.toml");
        std::fs::write(&cfg_path, "").unwrap();

        let cfg = SimConfig {
            network: Some("nets/tiny.toml".into()),
            ..SimConfig::default()
        };
        let spec = cfg
            .load_network(cfg_path.to_str())
            .expect("resolved next to the configuration");

        assert_eq!(spec.nodes.len(), 1);
        let _ = std::fs::remove_dir_all(&dir);
    }

    #[test]
    fn a_missing_network_file_is_reported_with_its_path() {
        let cfg = SimConfig {
            network: Some("nowhere/at/all.toml".into()),
            ..SimConfig::default()
        };

        let Err(CellweaveError::Config(message)) = cfg.load_network(None) else {
            panic!("a missing network must be refused");
        };
        assert!(message.contains("nowhere/at/all.toml"), "{message}");
    }
}
