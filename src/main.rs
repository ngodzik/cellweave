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
use cellweave_core::{
    CellId, CellInputs, CellKind, Concentration, DiffusionField, Pos2, SimSnapshot, TimeStep,
};
use cellweave_engine::coupling::{apply_uptake, hold_medium_bath, mean_per_cell};
use cellweave_engine::{CpmLattice, RasErkNetwork, ScalarField};
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

/// Oxygen diffusion coefficient, in pixels squared per unit time.
const OXYGEN_DIFFUSION: f64 = 0.2;

/// What one living pixel takes up per unit time, as a fraction of the bath.
///
/// With the bath at 1.0 this sets how deep a tissue can be before its middle
/// starves: a bathed disc of radius R sits q R² / (4 D) below the surface at
/// its centre, so with D = 0.2 the centre reaches the death level near R = 50
/// pixels, about thirty cells. Not calibrated. A real spheroid keeps a viable
/// rim several cells thick; this value gives about one, and a lower one would
/// push the core past the edge of the default grid before it forms.
const OXYGEN_UPTAKE: f64 = 0.0003;

/// The field is relaxed until no pixel moves by more than this per step.
const OXYGEN_TOLERANCE: f64 = 1e-4;

/// Give up on relaxing the field after this many steps, loudly.
const OXYGEN_MAX_STEPS: usize = 20_000;

/// Units of network time that pass in one lattice step, taken in sub steps of
/// one unit each so the RK4 scheme stays well inside its stability limit.
///
/// The three clocks have to agree. The field is relaxed to equilibrium every
/// lattice step, so it is instantaneous. A cell grows to division in about
/// twenty lattice steps. The network settles in about fifteen units of its own
/// time, so at one unit per lattice step it would answer barely faster than the
/// cell cycle, and a cell's fate would trail its surroundings by most of a
/// generation. Eight units per step brings the answer down to two lattice
/// steps: field, then network, then cycle, each well ahead of the next, which is
/// the order biology has them in even if the ratios here are compressed.
const NETWORK_TIME_PER_MCS: u32 = 8;

/// Survival signal below which a cell is necrotic.
const DEATH_THRESHOLD: f64 = 0.25;

/// Hypoxia response above which a cell stops dividing while staying alive.
const QUIESCENCE_THRESHOLD: f64 = 0.45;

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

    // Oxygen, full everywhere to start with, bathed wherever no cell sits.
    let mut oxygen = ScalarField::new(config.width, config.height, OXYGEN_DIFFUSION, 0.0, 1.0)
        .map_err(|e| anyhow::anyhow!("{e}"))?;

    let mut rng = StdRng::seed_from_u64(42);

    for mcs in 0..config.mcs {
        // 1. Downward coupling: the field settles against the current layout of
        //    the cells, and each cell reads the oxygen where it actually sits.
        apply_uptake(&mut oxygen, &lattice, OXYGEN_UPTAKE);
        hold_medium_bath(&mut oxygen, &lattice, 1.0);
        let relax_steps = oxygen
            .relax(TimeStep(1.0), OXYGEN_TOLERANCE, OXYGEN_MAX_STEPS)
            .map_err(|e| anyhow::anyhow!("at MCS {mcs}: {e}"))?;
        let local_o2 = mean_per_cell(&oxygen, &lattice);

        // 2. Advance every living network by one unit of time. A dead cell's
        //    network is frozen where it stopped.
        for (i, net) in networks.iter_mut().enumerate() {
            if lattice.cell_kind(cell_of(i)) == Some(CellKind::Necrotic) {
                continue;
            }
            let inputs = CellInputs {
                o2: Concentration(local_o2[i + 1]),
                egf: Concentration(0.2),
                vegf: Concentration(0.0),
                drugs: vec![],
            };
            for _ in 0..NETWORK_TIME_PER_MCS {
                net.step(TimeStep(1.0), &inputs);
            }
        }

        // 3. Fate and upward coupling: a cell whose survival collapsed is
        //    necrotic and holds its shape; every other cell's network decides
        //    how large it tries to be.
        for (i, net) in networks.iter().enumerate() {
            let id = cell_of(i);
            match lattice.cell_kind(id) {
                Some(CellKind::Necrotic) => {}
                Some(_) if net.survival() < DEATH_THRESHOLD => {
                    lattice.set_kind(id, CellKind::Necrotic);
                    if let Some(v) = lattice.volume(id) {
                        lattice.set_volume_target(id, v, net.mechanics().lambda_volume);
                    }
                }
                Some(_) => {
                    let m = net.mechanics();
                    lattice.set_volume_target(id, m.target_volume, m.lambda_volume);
                }
                None => {}
            }
        }

        // 4. Advance CPM.
        lattice.monte_carlo_step(&mut rng);

        // The bath is defined by what it can reach from the edge of the box. The
        // moment tissue touches that edge, the box is deciding what the tissue
        // sees, and everything after this point would be an artefact of its
        // size. Stop, and say so.
        if touches_the_edge(&lattice) {
            eprintln!(
                "stopped early: tissue reached the edge of the box at MCS {mcs} with {} cells. \
                 Results up to here are valid; to go further, use a larger grid.",
                networks.len()
            );
            break;
        }

        // 5. Cells that have grown enough, and are not too short of oxygen to
        //    cycle, split in two.
        if mcs % DIVISION_CHECK_EVERY == 0 && networks.len() < MAX_CELLS {
            divide_ready_cells(&mut lattice, &mut networks, &mut rng);
        }

        // 6. Save snapshot at requested interval. Daughters born this step have
        //    not read the field yet; dividing only relabels pixels, so reading it
        //    again against the new layout is exact.
        if mcs % config.snapshot_interval == 0 {
            let local_o2 = mean_per_cell(&oxygen, &lattice);
            let cells = lattice
                .cell_records()
                .zip(networks.iter())
                .map(|((id, rec), net)| cellweave_core::CellSnapshot {
                    id,
                    kind: rec.kind,
                    volume: rec.volume,
                    proteins: net.state().clone(),
                    oxygen: Concentration(local_o2[id.0 as usize]),
                })
                .collect();

            let snapshot = SimSnapshot { mcs, cells };
            output
                .write_snapshot(&snapshot)
                .map_err(|e| anyhow::anyhow!("{e}"))?;

            let necrotic = lattice
                .cell_records()
                .filter(|(_, r)| r.kind == CellKind::Necrotic)
                .count();
            eprintln!(
                "  MCS {mcs:5} | cells {:4} | necrotic {necrotic:4} | centre O2 {:.3} | field settled in {relax_steps:4} steps",
                networks.len(),
                oxygen.concentration_at(cx, cy)
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

/// Split every tumor cell that has grown past the division threshold and is
/// still cycling.
///
/// Hypoxia arrests the cell cycle long before it kills: a cell whose hypoxia
/// response is high is quiescent, alive and holding its place, but not
/// dividing. Without that middle state nothing bounds growth short of death.
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
            let cycling = net.hypoxia_response() < QUIESCENCE_THRESHOLD;
            big_enough && cycling && lattice.cell_kind(id) == Some(CellKind::Tumor)
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

/// Whether any cell occupies a pixel on the edge of the box.
fn touches_the_edge(lattice: &CpmLattice) -> bool {
    let (w, h) = lattice.dims();
    (0..w).any(|x| lattice.occupant(x, 0) != 0 || lattice.occupant(x, h - 1) != 0)
        || (0..h).any(|y| lattice.occupant(0, y) != 0 || lattice.occupant(w - 1, y) != 0)
}
