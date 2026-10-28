//! cellweave: multiscale cancer cell simulator.
//!
//! Usage:
//!
//! ```text
//! cellweave [--config <path>] [--mcs <n>] [--seed <n>] [--output <dir>]
//!           [--suppress <clone>,<clone>...] [--attribute]
//!           [--seeds <a>-<b> [--versus <path>]]
//! cellweave --version
//! ```
//!
//! The loop itself is [`cellweave_engine::Simulation`]; this binary reads a
//! configuration, runs it, and writes what it sees.

#![deny(unsafe_code)]

mod attribute;
mod campaign;

use anyhow::{Context, Result};
use cellweave_core::Seeding;
use cellweave_core::traits::SimOutput;
use cellweave_engine::{Measurements, Simulation};
use cellweave_io::{JsonOutput, MeasurementLog, SimConfig, write_clones};

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
    // --seed and --output, so that one configuration can be run under many
    // seeds side by side without editing it.
    if let Some(w) = args.windows(2).find(|w| w[0] == "--seed") {
        config.seed = w[1].parse().context("--seed must be an integer")?;
    }
    if let Some(w) = args.windows(2).find(|w| w[0] == "--output") {
        config.output_dir = w[1].clone();
    }
    // --suppress: clones whose founding mutation this run draws but does not
    // make, to see what would have happened without them.
    if let Some(w) = args.windows(2).find(|w| w[0] == "--suppress") {
        config.suppress = w[1]
            .split(',')
            .map(|c| c.trim().parse())
            .collect::<Result<_, _>>()
            .context("--suppress takes clone numbers separated by commas")?;
    }

    // --seeds: the configuration under every seed of a range, and with
    // --versus, a second configuration under the same seeds, compared pair by
    // pair.
    if let Some(w) = args.windows(2).find(|w| w[0] == "--seeds") {
        let seeds = campaign::parse_seeds(&w[1])?;
        let arm = |config: &SimConfig, path: Option<&str>| -> Result<campaign::Arm> {
            Ok(campaign::Arm {
                name: path
                    .and_then(|p| std::path::Path::new(p).file_stem())
                    .map_or_else(
                        || "default".to_owned(),
                        |s| s.to_string_lossy().into_owned(),
                    ),
                scenario: config.scenario(),
                spec: config.load_network(path)?,
                steps: config.mcs,
            })
        };
        let mut arms = vec![arm(&config, config_path)?];
        if let Some(v) = args.windows(2).find(|w| w[0] == "--versus") {
            let mut other =
                SimConfig::from_file(&v[1]).with_context(|| format!("loading config {}", v[1]))?;
            other.mcs = config.mcs;
            arms.push(arm(&other, Some(&v[1]))?);
        }
        eprintln!(
            "campaign: {} arm(s), seeds {}..={}, {} MCS each, all at once",
            arms.len(),
            seeds.start(),
            seeds.end(),
            config.mcs
        );
        let summary = campaign::campaign(&arms, seeds, std::path::Path::new(&config.output_dir))?;
        eprint!("{summary}");
        return Ok(());
    }

    if args.iter().any(|a| a == "--attribute") {
        let spec = config.load_network(config_path)?;
        eprintln!(
            "attributing the outcome of {} MCS to the mutations of its leading clone, \
             one run per mutation",
            config.mcs
        );
        let table = attribute::attribute(
            &config.scenario(),
            &spec,
            config.mcs,
            std::path::Path::new(&config.output_dir),
        )?;
        eprint!("{table}");
        return Ok(());
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
                "  MCS {:5} | living {:4} | necrotic {:4} | born {:5} | gone {:5} | O2 at first {:.3} | settled {:4}{}",
                step.mcs,
                sim.living(),
                sim.necrotic(),
                sim.born(),
                sim.gone(),
                step.first_cell_oxygen,
                step.relax_steps,
                if step.exposure > 0.0 {
                    format!(" | stress {:.2}", step.exposure)
                } else {
                    String::new()
                }
            );
        }
    }

    output.finalize()?;
    log.finish()?;
    let sizes = sim.clone_sizes();
    write_clones(&config.output_dir, sim.clones(), &sizes)?;
    if let Some(summary) = summarise(&measured) {
        eprintln!("{summary}");
    }
    let steps = measured.len() as u64;
    if let Some(days) = config.days(steps)? {
        eprintln!(
            "  {steps} steps at {} hours each: about {days:.0} days",
            config.hours_per_step.unwrap_or_default()
        );
    }
    if let Some(sisters) = sisters(&sim) {
        eprintln!("{sisters}");
    }
    if matches!(sim.scenario().seeding, Seeding::Sheet { .. })
        && let Some(scaling) = scaling(&sizes)
    {
        eprintln!("{scaling}");
    }
    if sim.clones().len() > 1 {
        eprintln!("{}", leading_clones(&sim, &sizes));
    }
    eprintln!("done.");
    Ok(())
}

/// Mutations shown per clone in the summary, the latest ones.
const SHOWN_MUTATIONS: usize = 4;

/// The clones holding the most living cells, each with the chain of
/// mutations that separates it from the founder.
fn leading_clones<S: cellweave_engine::simulation::SimState>(
    sim: &Simulation<S>,
    sizes: &[usize],
) -> String {
    let living: usize = sizes.iter().sum();
    let mut order: Vec<usize> = (0..sizes.len()).filter(|&c| sizes[c] > 0).collect();
    order.sort_by_key(|&c| std::cmp::Reverse(sizes[c]));
    let mut text = format!(
        "{} clones founded, {} still living",
        sim.clones().len(),
        order.len()
    );
    let mut carried: Vec<(String, f64)> = sim
        .coverage()
        .into_iter()
        .filter(|(_, share)| *share > 0.0)
        .collect();
    carried.sort_by(|a, b| b.1.total_cmp(&a.1));
    if !carried.is_empty() {
        let shares: Vec<String> = carried
            .iter()
            .map(|(node, share)| format!("{node} {:.0}%", 100.0 * share))
            .collect();
        text.push_str(&format!(
            "\n  the living carrying a change to: {}",
            shares.join(", ")
        ));
    }
    let sheet = matches!(sim.scenario().seeding, Seeding::Sheet { .. });
    for &c in order.iter().take(5) {
        let mut chain = Vec::new();
        let mut at = sim.clones().get(c);
        while let Some(clone) = at {
            if let Some(m) = &clone.mutation {
                chain.push(m.to_string());
            }
            at = clone.parent.and_then(|p| sim.clones().get(p.0 as usize));
        }
        chain.reverse();
        // A long line of descent is shown by its latest changes, the ones
        // that set it apart from its nearest relatives.
        let earlier = chain.len().saturating_sub(SHOWN_MUTATIONS);
        let latest = chain[earlier..].join(", then ");
        let chain = match (chain.is_empty(), earlier) {
            // In a sheet every cell founded a clone of its own.
            (true, _) if sheet => "a founder, unmutated".to_owned(),
            (true, _) => "the founder".to_owned(),
            (false, 0) => latest,
            (false, n) => format!("{n} earlier, then {latest}"),
        };
        text.push_str(&format!(
            "\n  clone {c:4}: {:3} cells ({:3.0}%), {chain}",
            sizes[c],
            100.0 * sizes[c] as f64 / living.max(1) as f64
        ));
    }
    text
}

/// What the run came to, in a few lines: the shape the tissue ended in, and
/// how fast it was turning over in the second half of the run, which leaves
/// out the growth from a single cell.
fn summarise(measured: &[Measurements]) -> Option<String> {
    let last = measured.last()?;
    let first = &measured[measured.len() / 2];
    let span = (last.mcs - first.mcs).max(1) as f64;
    let per_hundred = |n: u64| 100.0 * n as f64 / span;
    // How often a living cell divides, on average: the living, averaged over
    // the window, times its length, over the births in it.
    let window = &measured[measured.len() / 2..];
    let living = window
        .iter()
        .map(|m| (m.cycling + m.quiescent) as f64)
        .sum::<f64>()
        / window.len() as f64;
    let births = last.born - first.born;
    let interval = if births > 0 {
        format!(
            "\n  a living cell divided once every {:.0} steps, on average",
            living * span / births as f64
        )
    } else {
        String::new()
    };
    Some(format!(
        "at MCS {}: {} cycling, {} quiescent, {} necrotic, {} apoptotic\n  \
         tissue radius {:.1} px, necrotic radius {:.1} px, viable depth {:.1} px\n  \
         from MCS {}: {:.1} births and {:.1} cells gone per 100 MCS{}",
        last.mcs,
        last.cycling,
        last.quiescent,
        last.necrotic,
        last.apoptotic,
        last.tissue_radius,
        last.necrotic_radius,
        last.viable_depth,
        first.mcs,
        per_hundred(last.born - first.born),
        per_hundred(last.gone - first.gone),
        interval,
    ))
}

/// What pairs of sisters did next, and the `r` of the single progenitor
/// model it gives, or `None` before a pair is complete.
fn sisters<S: cellweave_engine::simulation::SimState>(sim: &Simulation<S>) -> Option<String> {
    let fates = sim.sister_fates();
    let r = fates.r()?;
    let share = |n: u64| 100.0 * n as f64 / fates.pairs() as f64;
    Some(format!(
        "sisters, {} pairs: both divided {:.1}%, one of each {:.1}%, both left {:.1}%; \
         r = {r:.3}, a quarter if sisters were independent",
        fates.pairs(),
        share(fates.both_divided),
        share(fates.one_each),
        share(fates.both_left),
    ))
}

/// Multiples of the mean size the surviving clones are compared against.
const SCALING_MULTIPLES: [f64; 3] = [1.0, 2.0, 3.0];

/// How the surviving clones' sizes spread around their mean: the share
/// larger than one, two and three times it, beside the exponential that
/// neutral competition in a renewing epithelium settles into, or `None`
/// with no clone alive.
fn scaling(sizes: &[usize]) -> Option<String> {
    let surviving: Vec<f64> = sizes
        .iter()
        .filter(|&&s| s > 0)
        .map(|&s| s as f64)
        .collect();
    if surviving.is_empty() {
        return None;
    }
    let mean = surviving.iter().sum::<f64>() / surviving.len() as f64;
    let shares: Vec<String> = SCALING_MULTIPLES
        .iter()
        .map(|&x| {
            let larger = surviving.iter().filter(|&&s| s > x * mean).count();
            format!(
                "{:.2} (exponential {:.2})",
                larger as f64 / surviving.len() as f64,
                (-x).exp()
            )
        })
        .collect();
    Some(format!(
        "{} surviving clones, {:.1} cells each on average; larger than 1, 2, 3 times that: {}",
        surviving.len(),
        mean,
        shares.join(", ")
    ))
}
