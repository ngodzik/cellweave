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
use cellweave_engine::Simulation;
use cellweave_io::{JsonOutput, SimConfig};

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

    for _ in 0..config.mcs {
        let step = sim.step()?;

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
    eprintln!("done.");
    Ok(())
}
