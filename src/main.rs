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
    CellId, CellInputs, CellKind, Concentration, Pos2, SimSnapshot, TimeStep, Volume,
};
use cellweave_engine::coupling::{
    apply_background_uptake, apply_uptake, hold_medium_bath, hold_obstacles_at, mean_per_cell,
};
use cellweave_engine::{CpmLattice, Network, ScalarField};
use cellweave_io::{JsonOutput, OxygenSource, SimConfig};
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

/// The field is relaxed until no pixel moves by more than this fraction of the
/// uptake per step.
///
/// Relative on purpose: the per step change is what diffusion has not yet
/// balanced against the uptake, so a tolerance above the uptake itself would
/// declare the field settled without having balanced anything.
const OXYGEN_TOLERANCE: f64 = 0.01;

/// Give up on relaxing the field after this many steps, loudly.
const OXYGEN_MAX_STEPS: usize = 20_000;

/// Before the first step the field is solved from scratch, and a per-step
/// tolerance is a poor guard against a slow front creeping out from a small
/// source: each step moves less than the tolerance while the total drifts for
/// hundreds of steps. Once, at the start, the field is settled far more
/// tightly so that every later step is an increment on a true equilibrium.
const OXYGEN_WARMUP_TOLERANCE: f64 = 1e-8;

/// And it is allowed to take its time doing so.
const OXYGEN_WARMUP_MAX_STEPS: usize = 2_000_000;

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

/// Pixels a necrotic cell gives up per lattice step until nothing is left.
///
/// Dead tissue is cleared, slowly, and its place taken by whatever grows into
/// it. Without this a tumour is bounded by its box; with it, by turnover.
const RESORPTION_RATE: u32 = 2;

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

    run(config, config_path)
}

fn run(config: SimConfig, config_path: Option<&str>) -> Result<()> {
    eprintln!(
        "cellweave | grid {}×{} | {} MCS | {:?} | output → {}",
        config.width, config.height, config.mcs, config.oxygen_source, config.output_dir
    );

    let mut output = JsonOutput::new(&config.output_dir).context("creating output directory")?;

    let mut lattice = CpmLattice::new(config.width, config.height, config.temperature)
        .map_err(|e| anyhow::anyhow!("{e}"))?;
    let cx = config.width / 2;
    let cy = config.height / 2;

    // Where the oxygen comes from decides where the first cell goes: at the
    // centre of a bath, or against the wall of a vessel.
    let first = match config.oxygen_source {
        OxygenSource::Bath => lattice.add_cell(CellKind::Tumor, Pos2 { x: cx, y: cy }, 8),
        OxygenSource::Vessel { radius } => {
            lattice.add_obstacle_disc(Pos2 { x: cx, y: cy }, radius);
            lattice.add_cell(
                CellKind::Tumor,
                Pos2 {
                    x: cx + radius + 8,
                    y: cy,
                },
                8,
            )
        }
    };
    let newborn_volume = lattice.volume(first).map_or(0, |v| v.0);

    // The hypothesis under test, read from a file rather than compiled in.
    let spec = config
        .load_network(config_path)
        .map_err(|e| anyhow::anyhow!("{e}"))?;
    let seed_network =
        Network::from_spec(&spec, newborn_volume).map_err(|e| anyhow::anyhow!("{e}"))?;
    if let Some(roles) = &spec.roles {
        eprintln!(
            "network | {} nodes | growth {} | survival {} | hypoxia {}",
            spec.nodes.len(),
            roles.growth,
            roles.survival,
            roles.hypoxia_response
        );
    }

    // One network per cell, or None once the cell is gone. Identifier `i + 1`
    // names the cell whose slot is `i`, which holds as long as a daughter slot is
    // pushed exactly when a daughter cell is created.
    let mut networks: Vec<Option<Network>> = vec![Some(seed_network)];

    let initial = match config.oxygen_source {
        OxygenSource::Bath => 1.0,
        OxygenSource::Vessel { .. } => 0.0,
    };
    let mut oxygen = ScalarField::new(config.width, config.height, OXYGEN_DIFFUSION, 0.0, initial)
        .map_err(|e| anyhow::anyhow!("{e}"))?;

    let mut rng = StdRng::seed_from_u64(42);
    let mut born = 1u64;
    let mut gone = 0u64;

    // Solve the field from scratch once, tightly, so the loop below only ever
    // has to nudge an equilibrium.
    declare_sources(&mut oxygen, &lattice, &config);
    let warmup = oxygen
        .relax(
            TimeStep(1.0),
            OXYGEN_WARMUP_TOLERANCE,
            OXYGEN_WARMUP_MAX_STEPS,
        )
        .map_err(|e| anyhow::anyhow!("settling the initial field: {e}"))?;
    eprintln!("  initial field settled in {warmup} steps");

    for mcs in 0..config.mcs {
        // 1. Downward coupling: the field settles against the current layout of
        //    the cells, and each cell reads the oxygen where it actually sits.
        declare_sources(&mut oxygen, &lattice, &config);
        let relax_steps = oxygen
            .relax(
                TimeStep(1.0),
                OXYGEN_TOLERANCE * config.oxygen_uptake,
                OXYGEN_MAX_STEPS,
            )
            .map_err(|e| anyhow::anyhow!("at MCS {mcs}: {e}"))?;
        let local_o2 = mean_per_cell(&oxygen, &lattice);

        // 2. Advance every living network by NETWORK_TIME_PER_MCS units of its
        //    time. A dead cell's network is frozen where it stopped.
        for (i, slot) in networks.iter_mut().enumerate() {
            let Some(net) = slot else { continue };
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

        // 3. Fate and upward coupling. A necrotic cell gives up pixels until
        //    none is left, and is then gone for good; a cell whose survival just
        //    collapsed becomes necrotic; every other cell's network decides how
        //    large it tries to be.
        for (i, slot) in networks.iter_mut().enumerate() {
            let id = cell_of(i);
            let Some(net) = slot else { continue };
            match lattice.cell_kind(id) {
                Some(CellKind::Necrotic) => {
                    let volume = lattice.volume(id).map_or(0, |v| v.0);
                    if volume == 0 {
                        *slot = None;
                        gone += 1;
                    } else {
                        let target = Volume(volume.saturating_sub(RESORPTION_RATE));
                        lattice.set_volume_target(id, target, net.mechanics().lambda_volume);
                    }
                }
                Some(_) if net.survival() < DEATH_THRESHOLD => {
                    lattice.set_kind(id, CellKind::Necrotic);
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

        // The bath is defined by what it can reach from the edge of the box, and
        // a cord should never get near it. The moment tissue touches that edge,
        // the box is deciding what the tissue sees, and everything after this
        // point would be an artefact of its size. Stop, and say so.
        if touches_the_edge(&lattice) {
            write_snapshot(&mut output, mcs, &lattice, &networks, &oxygen)?;
            eprintln!(
                "stopped early: tissue reached the edge of the box at MCS {mcs} with {} living cells. \
                 Results up to here are valid; to go further, use a larger grid.",
                living(&lattice)
            );
            break;
        }

        // 5. Cells that have grown enough, and are not too short of oxygen to
        //    cycle, split in two.
        if mcs % DIVISION_CHECK_EVERY == 0 && networks.len() < MAX_CELLS {
            born += divide_ready_cells(&mut lattice, &mut networks, &mut rng);
        }

        // 6. Save snapshot at requested interval.
        if mcs % config.snapshot_interval == 0 {
            write_snapshot(&mut output, mcs, &lattice, &networks, &oxygen)?;

            let necrotic = lattice
                .cell_records()
                .filter(|(_, r)| r.kind == CellKind::Necrotic && r.volume.0 > 0)
                .count();
            eprintln!(
                "  MCS {mcs:5} | living {:4} | necrotic {necrotic:4} | born {born:5} | gone {gone:5} | O2 at first {:.3} | settled {relax_steps:4}",
                living(&lattice),
                local_o2.get(first.0 as usize).copied().unwrap_or(0.0)
            );
        }
    }

    output.finalize().map_err(|e| anyhow::anyhow!("{e}"))?;
    eprintln!("done.");
    Ok(())
}

/// Tell the field who consumes and who supplies, from the current layout.
///
/// Around a bath, living cells consume and the medium the bath can reach is
/// held full. Around a vessel there is no empty liquid but stroma, so the
/// medium consumes too, and the vessel wall is what is held.
fn declare_sources(oxygen: &mut ScalarField, lattice: &CpmLattice, config: &SimConfig) {
    apply_uptake(oxygen, lattice, config.oxygen_uptake);
    oxygen.release_all();
    match config.oxygen_source {
        OxygenSource::Bath => hold_medium_bath(oxygen, lattice, 1.0),
        OxygenSource::Vessel { .. } => {
            apply_background_uptake(oxygen, lattice, config.oxygen_uptake);
            hold_obstacles_at(oxygen, lattice, 1.0);
        }
    }
}

/// The cell whose slot is `i`. Identifier 0 is the medium.
fn cell_of(i: usize) -> CellId {
    CellId(i as u32 + 1)
}

/// Cells that are alive: not necrotic, and still holding pixels.
fn living(lattice: &CpmLattice) -> usize {
    lattice
        .cell_records()
        .filter(|(_, r)| r.kind != CellKind::Necrotic && r.volume.0 > 0)
        .count()
}

/// Split every tumor cell that has grown past the division threshold and is
/// still cycling, and say how many daughters were born.
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
    networks: &mut Vec<Option<Network>>,
    rng: &mut StdRng,
) -> u64 {
    let ready: Vec<usize> = networks
        .iter()
        .enumerate()
        .filter_map(|(i, slot)| slot.as_ref().map(|net| (i, net)))
        .filter(|(i, net)| {
            let id = cell_of(*i);
            let threshold = DIVISION_RATIO * f64::from(net.reference_volume());
            let big_enough = lattice
                .volume(id)
                .is_some_and(|v| f64::from(v.0) >= threshold);
            let cycling = net.hypoxia_response() < QUIESCENCE_THRESHOLD;
            big_enough && cycling && lattice.cell_kind(id) == Some(CellKind::Tumor)
        })
        .map(|(i, _)| i)
        .collect();

    let mut born = 0;
    for i in ready {
        if networks.len() >= MAX_CELLS {
            break;
        }
        let angle = rng.random::<f64>() * std::f64::consts::PI;
        if lattice.divide(cell_of(i), angle).is_some() {
            let daughter = networks[i].clone();
            networks.push(daughter);
            born += 1;
        }
    }

    // Kept as a hard assertion rather than a comment: were the two to drift, a
    // network would quietly drive the wrong cell and every result after that
    // would be wrong without a single test failing.
    assert_eq!(
        networks.len(),
        lattice.cell_count(),
        "one slot per cell: the lattice and the network list have diverged"
    );
    born
}

/// Whether any cell occupies a pixel on the edge of the box.
fn touches_the_edge(lattice: &CpmLattice) -> bool {
    let (w, h) = lattice.dims();
    let is_cell = |x: u32, y: u32| {
        let id = lattice.occupant(x, y);
        id != 0 && id != cellweave_engine::cpm::OBSTACLE
    };
    (0..w).any(|x| is_cell(x, 0) || is_cell(x, h - 1))
        || (0..h).any(|y| is_cell(0, y) || is_cell(w - 1, y))
}

/// Record every cell still holding pixels, with the oxygen it sits in.
///
/// Daughters born this step have not read the field yet; dividing only relabels
/// pixels, so reading the field again against the current layout is exact. A
/// cell that is gone has nothing left to record.
fn write_snapshot(
    output: &mut JsonOutput,
    mcs: u64,
    lattice: &CpmLattice,
    networks: &[Option<Network>],
    oxygen: &ScalarField,
) -> Result<()> {
    let local_o2 = mean_per_cell(oxygen, lattice);
    let cells = lattice
        .cell_records()
        .zip(networks.iter())
        .filter_map(|((id, rec), slot)| {
            let net = slot.as_ref()?;
            (rec.volume.0 > 0).then(|| cellweave_core::CellSnapshot {
                id,
                kind: rec.kind,
                volume: rec.volume,
                proteins: net.state().clone(),
                oxygen: Concentration(local_o2[id.0 as usize]),
            })
        })
        .collect();
    output
        .write_snapshot(&SimSnapshot { mcs, cells })
        .map_err(|e| anyhow::anyhow!("{e}"))
}
