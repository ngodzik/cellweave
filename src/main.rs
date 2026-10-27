//! cellweave: multiscale cancer cell simulator.
//!
//! Usage:
//!
//! ```text
//! cellweave [--config <path>] [--mcs <n>]
//! cellweave --version
//! ```
//!
//! The loop itself is [`cellweave_engine::Simulation`]; this binary reads a
//! configuration, runs it, and writes what it sees.

#![deny(unsafe_code)]

use anyhow::{Context, Result};
use cellweave_core::traits::SimOutput;
use cellweave_engine::{Measurements, Simulation};
use cellweave_io::{JsonOutput, MeasurementLog, SimConfig};

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
    let mut log = MeasurementLog::new(&config.output_dir)?;

    // The hypothesis under test, read from a file rather than compiled in.
    let spec = config.load_network(config_path)?;
    if let Some(roles) = &spec.roles {
        eprintln!(
            "network | {} nodes | growth {} | survival {} | hypoxia {}",
            spec.nodes.len(),
            roles.growth,
            roles.survival,
            roles.hypoxia_response
        );
    }

    let (mut sim, warmup) = Simulation::new(config.scenario(), &spec)?.start()?;
    eprintln!("  initial field settled in {warmup} steps");

    // Every step is measured, and the measurements kept, so that a rate can
    // be read off any window of the run once it is over.
    let mut measured = Vec::new();
    for _ in 0..config.mcs {
        let step = sim.step()?;
        let m = sim.measure();
        log.record(&m)?;
        measured.push(m);

        // The moment tissue touches the edge of the box, the box is deciding
        // what the tissue sees. Stop, and say so.
        if step.reached_edge {
            output.write_snapshot(&sim.snapshot())?;
            eprintln!(
                "stopped early: tissue reached the edge of the box at MCS {} with {} living cells. \
                 Results up to here are valid; to go further, use a larger grid.",
                step.mcs,
                sim.living()
            );
            break;
        }

        if step.mcs.is_multiple_of(config.snapshot_interval) {
            output.write_snapshot(&sim.snapshot())?;
            eprintln!(
                "  MCS {:5} | living {:4} | necrotic {:4} | born {:5} | gone {:5} | O2 at first {:.3} | settled {:4}",
                step.mcs,
                sim.living(),
                sim.necrotic(),
                sim.born(),
                sim.gone(),
                step.first_cell_oxygen,
                step.relax_steps
            );
        }
    }

    output.finalize()?;
    log.finish()?;
    if let Some(summary) = summarise(&measured) {
        eprintln!("{summary}");
    }
    eprintln!("done.");
    Ok(())
}

/// What the run came to, in a few lines: the shape the tissue ended in, and
/// how fast it was turning over in the second half of the run, which leaves
/// out the growth from a single cell.
fn summarise(measured: &[Measurements]) -> Option<String> {
    let last = measured.last()?;
    let first = &measured[measured.len() / 2];
    let span = (last.mcs - first.mcs).max(1) as f64;
    let per_hundred = |n: u64| 100.0 * n as f64 / span;
    Some(format!(
        "at MCS {}: {} cycling, {} quiescent, {} necrotic\n  \
         tissue radius {:.1} px, necrotic radius {:.1} px, viable depth {:.1} px\n  \
         from MCS {}: {:.1} births and {:.1} resorptions per 100 MCS",
        last.mcs,
        last.cycling,
        last.quiescent,
        last.necrotic,
        last.tissue_radius,
        last.necrotic_radius,
        last.viable_depth,
        first.mcs,
        per_hundred(last.born - first.born),
        per_hundred(last.gone - first.gone),
    ))
}
