//! Simulation output: JSON snapshots, and a table of measurements, to disk.

use cellweave_core::{CellweaveError, SimOutput, SimSnapshot};
use cellweave_engine::{CloneRecord, Measurements};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};

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

/// Writes one row of [`Measurements`] per step into `measurements.csv`.
///
/// A table rather than JSON so that it opens as it is in a spreadsheet or a
/// plotting library, and one row per step rather than per snapshot because the
/// numbers are cheap and a turnover rate needs every step to be read off.
/// Floats are written to four decimals, a fixed width, so that two runs of
/// one seed write the same bytes.
pub struct MeasurementLog {
    out: BufWriter<File>,
}

/// The columns, in the order [`MeasurementLog::record`] writes them.
pub const MEASUREMENT_COLUMNS: &str = "mcs,cycling,quiescent,necrotic,differentiating,born,gone,clones,tissue_radius,necrotic_radius,viable_depth";

impl MeasurementLog {
    /// Start the table in `dir`, which must exist, with its header row.
    ///
    /// # Errors
    ///
    /// [`CellweaveError::Simulation`] when the file cannot be written.
    pub fn new(dir: impl AsRef<Path>) -> Result<Self, CellweaveError> {
        let path = dir.as_ref().join("measurements.csv");
        let file = File::create(&path)
            .map_err(|e| CellweaveError::Simulation(format!("{}: {e}", path.display())))?;
        let mut log = Self {
            out: BufWriter::new(file),
        };
        log.line(MEASUREMENT_COLUMNS)?;
        Ok(log)
    }

    /// Append one step.
    ///
    /// # Errors
    ///
    /// [`CellweaveError::Simulation`] when the file cannot be written.
    pub fn record(&mut self, m: &Measurements) -> Result<(), CellweaveError> {
        self.line(&format!(
            "{},{},{},{},{},{},{},{},{:.4},{:.4},{:.4}",
            m.mcs,
            m.cycling,
            m.quiescent,
            m.necrotic,
            m.differentiating,
            m.born,
            m.gone,
            m.clones,
            m.tissue_radius,
            m.necrotic_radius,
            m.viable_depth
        ))
    }

    /// Write what is buffered.
    ///
    /// # Errors
    ///
    /// [`CellweaveError::Simulation`] when the file cannot be written.
    pub fn finish(mut self) -> Result<(), CellweaveError> {
        self.out
            .flush()
            .map_err(|e| CellweaveError::Simulation(e.to_string()))
    }

    fn line(&mut self, text: &str) -> Result<(), CellweaveError> {
        writeln!(self.out, "{text}").map_err(|e| CellweaveError::Simulation(e.to_string()))
    }
}

/// Write every clone ever founded to `clones.json` in `dir`: its parent, the
/// step it arose at, the mutation that founded it, and how many living cells
/// it has in `sizes`, indexed like the clones.
///
/// # Errors
///
/// [`CellweaveError::Simulation`] when the file cannot be written.
pub fn write_clones(
    dir: impl AsRef<Path>,
    clones: &[CloneRecord],
    sizes: &[usize],
) -> Result<(), CellweaveError> {
    let list: Vec<serde_json::Value> = clones
        .iter()
        .map(|c| {
            serde_json::json!({
                "id": c.id.0,
                "parent": c.parent.map(|p| p.0),
                "founded_at": c.founded_at,
                "mutation": c.mutation.as_ref().map(|m| serde_json::json!({
                    "node": m.node,
                    "term": m.term,
                    "kind": m.kind,
                    "factor": m.factor,
                })),
                "living": sizes.get(c.id.0 as usize).copied().unwrap_or(0),
            })
        })
        .collect();
    let json = serde_json::to_string_pretty(&list)
        .map_err(|e| CellweaveError::Simulation(e.to_string()))?;
    fs::write(dir.as_ref().join("clones.json"), json)
        .map_err(|e| CellweaveError::Simulation(e.to_string()))
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]

    use super::*;

    #[test]
    fn the_table_has_a_header_and_one_fixed_width_row_per_step() {
        let dir = std::env::temp_dir().join(format!("cellweave-log-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let mut log = MeasurementLog::new(&dir).unwrap();
        log.record(&Measurements {
            mcs: 7,
            cycling: 3,
            quiescent: 2,
            necrotic: 1,
            differentiating: 0,
            born: 9,
            gone: 4,
            clones: 2,
            tissue_radius: 12.0,
            necrotic_radius: 1.0 / 3.0,
            viable_depth: 11.666_666,
        })
        .unwrap();
        log.finish().unwrap();

        let text = fs::read_to_string(dir.join("measurements.csv")).unwrap();
        fs::remove_dir_all(&dir).unwrap();
        let lines: Vec<&str> = text.lines().collect();
        assert_eq!(lines[0], MEASUREMENT_COLUMNS);
        assert_eq!(lines[1], "7,3,2,1,0,9,4,2,12.0000,0.3333,11.6667");
        assert_eq!(
            lines[0].split(',').count(),
            lines[1].split(',').count(),
            "every column has a value"
        );
    }
}
