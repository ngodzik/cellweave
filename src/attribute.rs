//! Which events made a history.
//!
//! Run the configuration once as it is, find the clone holding the most
//! living cells at the end, and walk back through the mutations that made it.
//! Then run the same configuration again once per mutation, with that one
//! mutation drawn but not made, everything else identical, and say what
//! changed. A mutation that changes nothing leaves the run identical to the
//! last bit, which no amount of noise can fake; one that does shows when the
//! two histories parted, and where the second one ended.

use std::fmt::Write as _;
use std::fs;
use std::path::Path;

use anyhow::{Context, Result};
use cellweave_core::{CloneId, NetworkSpec};
use cellweave_engine::simulation::Scenario;
use cellweave_engine::{CloneRecord, Measurements, Simulation};
use rayon::prelude::*;

/// What one run came to, and enough of its path to compare it with another.
struct History {
    /// The tissue measured at every step taken.
    measured: Vec<Measurements>,
    /// Who held each pixel at the end, row by row.
    lattice: Vec<u32>,
    /// Every clone founded, and its living cells at the end.
    clones: Vec<CloneRecord>,
    sizes: Vec<usize>,
    /// The run stopped at the edge of the box.
    reached_edge: bool,
}

/// Run `scenario` for at most `steps` steps, stopping at the edge.
fn run(scenario: Scenario, spec: &NetworkSpec, steps: u64) -> Result<History> {
    let (mut sim, _) = Simulation::new(scenario, spec)?.start()?;
    let mut measured = Vec::new();
    let mut reached_edge = false;
    for _ in 0..steps {
        let step = sim.step()?;
        measured.push(sim.measure());
        if step.reached_edge {
            reached_edge = true;
            break;
        }
    }
    let lattice = sim.lattice();
    let (w, h) = lattice.dims();
    let grid = (0..h)
        .flat_map(|y| (0..w).map(move |x| (x, y)))
        .map(|(x, y)| lattice.occupant(x, y))
        .collect();
    Ok(History {
        measured,
        lattice: grid,
        clones: sim.clones().to_vec(),
        sizes: sim.clone_sizes(),
        reached_edge,
    })
}

/// The mutations that separate `clone` from the founder, oldest first.
fn lineage(clones: &[CloneRecord], clone: CloneId) -> Vec<&CloneRecord> {
    let mut chain = Vec::new();
    let mut at = clones.get(clone.0 as usize);
    while let Some(record) = at {
        if record.mutation.is_some() {
            chain.push(record);
        }
        at = record.parent.and_then(|p| clones.get(p.0 as usize));
    }
    chain.reverse();
    chain
}

/// How a counterfactual history compares with the real one.
struct Comparison {
    clone: CloneId,
    founded_at: u64,
    mutation: String,
    /// The two histories are the same to the last bit: every measurement at
    /// every step, and every pixel at the end.
    identical: bool,
    /// The first step at which the measurements differ.
    parted_at: Option<u64>,
    /// Where the counterfactual history ended.
    end: Measurements,
    reached_edge: bool,
}

/// Attribute the outcome of the configured run to the mutations of its
/// leading clone, write `attribution.csv` into `dir`, and return the table as
/// text.
pub fn attribute(
    scenario: &Scenario,
    spec: &NetworkSpec,
    steps: u64,
    dir: &Path,
) -> Result<String> {
    let real = run(scenario.clone(), spec, steps)?;
    let Some(leader) = (0..real.sizes.len()).max_by_key(|&c| (real.sizes[c], std::cmp::Reverse(c)))
    else {
        return Ok("no clone is alive at the end; nothing to attribute".to_owned());
    };
    let leader = CloneId(leader as u32);
    let chain = lineage(&real.clones, leader);
    if chain.is_empty() {
        return Ok(format!(
            "the leading clone at the end, clone {}, carries no mutation; nothing to attribute",
            leader.0
        ));
    }

    let comparisons: Vec<Comparison> = chain
        .par_iter()
        .map(|record| -> Result<Comparison> {
            let mut without = scenario.clone();
            without.suppress.push(record.id);
            let other = run(without, spec, steps)?;
            let parted_at = real
                .measured
                .iter()
                .zip(&other.measured)
                .find(|(a, b)| a != b)
                .map(|(a, _)| a.mcs)
                .or_else(|| {
                    // Same as far as both went, but one went further.
                    (real.measured.len() != other.measured.len())
                        .then(|| real.measured.len().min(other.measured.len()) as u64)
                });
            let identical = parted_at.is_none() && real.lattice == other.lattice;
            let end = *other
                .measured
                .last()
                .context("a run takes at least one step")?;
            Ok(Comparison {
                clone: record.id,
                founded_at: record.founded_at,
                mutation: record
                    .mutation
                    .as_ref()
                    .map_or_else(String::new, ToString::to_string),
                identical,
                parted_at,
                end,
                reached_edge: other.reached_edge,
            })
        })
        .collect::<Result<_>>()?;

    let last = real
        .measured
        .last()
        .context("a run takes at least one step")?;
    let ending = |edge: bool, m: &Measurements| {
        format!(
            "{} at step {}: {} cycling, {} quiescent",
            if edge {
                "touched the edge"
            } else {
                "ran to the end"
            },
            m.mcs,
            m.cycling,
            m.quiescent
        )
    };
    let mut text = format!(
        "real history: leading clone {} ({} cells), {}\n",
        leader.0,
        real.sizes[leader.0 as usize],
        ending(real.reached_edge, last)
    );
    let mut csv = String::from(
        "clone,founded_at,mutation,identical,parted_at,reached_edge,end_mcs,cycling,quiescent,necrotic\n",
    );
    for c in &comparisons {
        let effect = if c.identical {
            "no effect: identical to the last bit".to_owned()
        } else if c.parted_at.is_none() {
            "the same measurements at every step, the cells laid out differently at the end"
                .to_owned()
        } else {
            format!(
                "parted at step {}; {}",
                c.parted_at.unwrap_or_default(),
                ending(c.reached_edge, &c.end)
            )
        };
        let _ = writeln!(
            text,
            "  without clone {:3} (step {:4}, {}): {effect}",
            c.clone.0, c.founded_at, c.mutation
        );
        let _ = writeln!(
            csv,
            "{},{},{},{},{},{},{},{},{},{}",
            c.clone.0,
            c.founded_at,
            c.mutation,
            c.identical,
            c.parted_at.map_or_else(String::new, |s| s.to_string()),
            c.reached_edge,
            c.end.mcs,
            c.end.cycling,
            c.end.quiescent,
            c.end.necrotic
        );
    }
    fs::create_dir_all(dir).with_context(|| format!("creating {}", dir.display()))?;
    fs::write(dir.join("attribution.csv"), csv).context("writing attribution.csv")?;
    Ok(text)
}
