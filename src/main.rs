//! cellweave: multiscale cancer cell simulator.
//!
//! Usage:
//!
//! ```text
//! cellweave [--config <path>] [--mcs <n>]
//! cellweave --version
//! ```

use anyhow::{Context, Result};
use cellweave_core::traits::SignalingNetwork;
use cellweave_core::traits::SimOutput;
use cellweave_core::{CellInputs, CellKind, Concentration, Pos2, SimSnapshot, TimeStep};
use cellweave_engine::{CpmLattice, RasErkNetwork};
use cellweave_io::{JsonOutput, SimConfig};
use rand::SeedableRng;
use rand::rngs::StdRng;

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();

    if args.iter().any(|a| a == "--version" || a == "-V") {
        println!("cellweave {}", env!("CARGO_PKG_VERSION"));
        return Ok(());
    }

    // Locate optional config file.
    let config_path = args
        .windows(2)
        .find(|w| w[0] == "--config")
        .map(|w| w[1].as_str());

    let mut config = match config_path {
        Some(p) => SimConfig::from_file(p).with_context(|| format!("loading config {p}"))?,
        None => SimConfig::default(),
    };

    // --mcs override.
    if let Some(w) = args.windows(2).find(|w| w[0] == "--mcs") {
        config.mcs = w[1].parse().context("--mcs must be an integer")?;
    }

    run(config)
}

fn run(config: SimConfig) -> Result<()> {
    eprintln!(
        "cellweave | grid {}×{} | {} MCS | output → {}",
        config.width, config.height, config.mcs, config.output_dir
    );

    let mut output = JsonOutput::new(&config.output_dir).context("creating output directory")?;

    // Build the CPM lattice.
    let mut lattice = CpmLattice::new(config.width, config.height, config.temperature)
        .map_err(|e| anyhow::anyhow!("{e}"))?;

    // Seed one tumor cell at the centre.
    let cx = config.width / 2;
    let cy = config.height / 2;
    let cell_id = lattice.add_cell(CellKind::Tumor, Pos2 { x: cx, y: cy }, 8);

    // One ODE network per cell.
    let mut networks: Vec<RasErkNetwork> = vec![RasErkNetwork::new(200)];

    let mut rng = StdRng::seed_from_u64(42);

    for mcs in 0..config.mcs {
        // 1. Advance ODE networks (1 second per MCS as a placeholder).
        let inputs = CellInputs {
            o2: Concentration(0.8),
            egf: Concentration(0.2),
            vegf: Concentration(0.0),
            drugs: vec![],
        };
        for net in networks.iter_mut() {
            net.step(TimeStep(1.0), &inputs);
        }

        // 2. Advance CPM.
        lattice.monte_carlo_step(&mut rng);

        // 3. Save snapshot at requested interval.
        if mcs % config.snapshot_interval == 0 {
            let cells = lattice
                .cell_records()
                .zip(networks.iter())
                .map(|((id, rec), net)| cellweave_core::CellSnapshot {
                    id,
                    kind: rec.kind,
                    volume: rec.volume,
                    proteins: net.state().clone(),
                })
                .collect();

            let snapshot = SimSnapshot { mcs, cells };
            output
                .write_snapshot(&snapshot)
                .map_err(|e| anyhow::anyhow!("{e}"))?;

            eprintln!(
                "  MCS {:5} | cell volume = {}",
                mcs,
                lattice.volume(cell_id).map_or(0, |v| v.0)
            );
        }
    }

    output.finalize().map_err(|e| anyhow::anyhow!("{e}"))?;
    eprintln!("done.");
    Ok(())
}
