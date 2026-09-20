//! cellweave: multiscale cancer cell simulator.
//!
//! Usage:
//!
//! ```text
//! cellweave [--config <path>] [--mcs <n>]
//! cellweave --version
//! ```

#![deny(unsafe_code)]

use anyhow::{Context, Result};
use cellweave_core::traits::SignalingNetwork;
use cellweave_core::traits::SimOutput;
use cellweave_core::{CellId, CellInputs, CellKind, Concentration, Pos2, SimSnapshot, TimeStep};
use cellweave_engine::{CpmLattice, RasErkNetwork};
use cellweave_io::{JsonOutput, SimConfig};
use rand::rngs::StdRng;
use rand::{RngExt, SeedableRng};

/// A cell divides once it reaches this multiple of its lineage's reference volume.
const DIVISION_RATIO: f64 = 1.4;

/// Monte Carlo steps between two checks for cells ready to divide.
///
/// Dividing scans the whole lattice, so it is not something to attempt on every
/// step, and a cell needs time to grow back anyway.
const DIVISION_CHECK_EVERY: u64 = 20;

/// Refuse to go past this many cells, as a guard against runaway growth.
const MAX_CELLS: usize = 2000;

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
    let first = lattice.add_cell(CellKind::Tumor, Pos2 { x: cx, y: cy }, 8);
    let newborn_volume = lattice.volume(first).map_or(0, |v| v.0);

    // One ODE network per cell. Identifier `i + 1` names the cell whose network
    // sits at index `i`, which holds as long as a daughter network is pushed
    // exactly when a daughter cell is created.
    let mut networks: Vec<RasErkNetwork> = vec![RasErkNetwork::new(newborn_volume)];

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

        // 2. Upward coupling: each network decides how large its cell tries to be.
        for (i, net) in networks.iter().enumerate() {
            let m = net.mechanics();
            lattice.set_volume_target(cell_of(i), m.target_volume, m.lambda_volume);
        }

        // 3. Advance CPM.
        lattice.monte_carlo_step(&mut rng);

        // 4. Cells that have grown enough split in two.
        if mcs % DIVISION_CHECK_EVERY == 0 && networks.len() < MAX_CELLS {
            divide_ready_cells(&mut lattice, &mut networks, &mut rng);
        }

        // 5. Save snapshot at requested interval.
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
                "  MCS {mcs:5} | cells {:4} | first cell volume {}",
                networks.len(),
                lattice.volume(first).map_or(0, |v| v.0)
            );
        }
    }

    output.finalize().map_err(|e| anyhow::anyhow!("{e}"))?;
    eprintln!("done.");
    Ok(())
}

/// The cell whose network sits at index `i`. Identifier 0 is the medium.
fn cell_of(i: usize) -> CellId {
    CellId(i as u32 + 1)
}

/// Split every living tumor cell that has grown past the division threshold.
///
/// A daughter inherits a copy of its parent's network, so its protein state and
/// its parameters, which is what makes the population a lineage rather than a
/// crowd of identical newborns.
fn divide_ready_cells(
    lattice: &mut CpmLattice,
    networks: &mut Vec<RasErkNetwork>,
    rng: &mut StdRng,
) {
    let ready: Vec<usize> = networks
        .iter()
        .enumerate()
        .filter(|(i, net)| {
            let id = cell_of(*i);
            let threshold = DIVISION_RATIO * f64::from(net.base_target_volume());
            let big_enough = lattice
                .volume(id)
                .is_some_and(|v| f64::from(v.0) >= threshold);
            big_enough && lattice.cell_kind(id) == Some(CellKind::Tumor)
        })
        .map(|(i, _)| i)
        .collect();

    for i in ready {
        if networks.len() >= MAX_CELLS {
            break;
        }
        let angle = rng.random::<f64>() * std::f64::consts::PI;
        if lattice.divide(cell_of(i), angle).is_some() {
            let daughter = networks[i].clone();
            networks.push(daughter);
        }
    }

    // Kept as a hard assertion rather than a comment: were the two to drift, a
    // network would quietly drive the wrong cell and every result after that
    // would be wrong without a single test failing.
    assert_eq!(
        networks.len(),
        lattice.cell_count(),
        "one network per cell: the lattice and the network list have diverged"
    );
}
