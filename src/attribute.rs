//! Which events made a history.
//!
//! Run the configuration once as it is, find the clone holding the most
//! living cells at the end, and walk back through the mutations that made it.
//! Then run the same configuration again once per mutation, with that one
//! mutation drawn but not made, everything else identical, and say what
//! changed. A mutation that changes nothing leaves the run identical to the
//! last bit, which no amount of noise can fake; one that does shows when the
//! two histories parted, and where the second one ended.
//!
//! When the configuration places a lesion, the question changes: not what
//! made the leading clone, but what the clones around the lesion did to it.
//! Each mutation carried by a clone that touched the lesion is undone at the
//! step the lesion appears, in every cell carrying it, and the history
//! replayed from there; until that step the two are the same run.

use std::collections::BTreeSet;
use std::fmt::Write as _;
use std::fs;
use std::path::Path;

use anyhow::{Context, Result};
use cellweave_core::{CellId, CloneId, NetworkSpec};
use cellweave_engine::simulation::{Reversion, Scenario, SimState};
use cellweave_engine::{CloneRecord, Fate, Measurements, Simulation};
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
    /// Clones with a living cell next to a living cell of the lesion, at
    /// any step while it lived, when there is one.
    touched: BTreeSet<CloneId>,
    /// Living cells in each clone at the end of the step the lesion
    /// appeared at.
    sizes_at_lesion: Vec<usize>,
    /// The first step at whose end no cell of the lesion was left.
    lesion_lost_at: Option<u64>,
}

/// Run `scenario` for at most `steps` steps, stopping at the edge.
fn run(scenario: Scenario, spec: &NetworkSpec, steps: u64) -> Result<History> {
    let lesion_at = scenario.lesion.as_ref().map(|l| l.at);
    let (mut sim, _) = Simulation::new(scenario, spec)?.start()?;
    let mut measured = Vec::new();
    let mut reached_edge = false;
    let mut touched = BTreeSet::new();
    let mut sizes_at_lesion = Vec::new();
    let mut lesion_lost_at = None;
    for _ in 0..steps {
        let step = sim.step()?;
        measured.push(sim.measure());
        if let Some(at) = lesion_at
            && step.mcs >= at
            && lesion_lost_at.is_none()
        {
            if step.mcs == at {
                sizes_at_lesion = sim.clone_sizes();
            }
            if sim.lesion_cells() == 0 {
                lesion_lost_at = Some(step.mcs);
            } else {
                touching(&sim, &mut touched);
            }
        }
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
        touched,
        sizes_at_lesion,
        lesion_lost_at,
    })
}

/// Add every clone with a living cell next to a living cell of the lesion.
fn touching<S: SimState>(sim: &Simulation<S>, touched: &mut BTreeSet<CloneId>) {
    let lattice = sim.lattice();
    // Per cell: `None` unless living; then whether it is the lesion's.
    let mut lesion: Vec<Option<bool>> = vec![None; lattice.cell_count() + 1];
    for (id, record) in lattice.cell_records() {
        if record.volume.0 > 0
            && matches!(sim.fate(id), Some(Fate::Cycling | Fate::Quiescent))
            && let Some(slot) = lesion.get_mut(id.0 as usize)
        {
            *slot = Some(sim.in_lesion(id));
        }
    }
    let (w, h) = lattice.dims();
    let kind = |id: u32| lesion.get(id as usize).copied().flatten();
    for y in 0..h {
        for x in 0..w {
            let a = lattice.occupant(x, y);
            let right = (x + 1 < w).then(|| lattice.occupant(x + 1, y));
            let below = (y + 1 < h).then(|| lattice.occupant(x, y + 1));
            for b in [right, below].into_iter().flatten() {
                if let (Some(ka), Some(kb)) = (kind(a), kind(b))
                    && ka != kb
                {
                    let other = if ka { b } else { a };
                    touched.extend(sim.clone_of(CellId(other)));
                }
            }
        }
    }
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
    if let Some(at) = scenario.lesion.as_ref().map(|l| l.at) {
        return attribute_lesion(scenario, spec, steps, dir, at);
    }
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

/// What a counterfactual did to the lesion, against the real history.
struct LesionOutcome {
    /// The clones undone, one or all of the candidates.
    undone: Vec<CloneId>,
    /// The two histories are the same to the last bit.
    identical: bool,
    /// The first step at which the measurements differ.
    parted_at: Option<u64>,
    /// Lesion cells at the end.
    lesion: usize,
    /// The step the lesion was lost at, if it was.
    lost_at: Option<u64>,
    /// The last step taken.
    end_mcs: u64,
}

/// Attribute what happened to the lesion to the mutations of the clones that
/// touched it, undoing each at the step the lesion appeared, then all of them
/// at once; write `attribution.csv` into `dir` and return the table as text.
fn attribute_lesion(
    scenario: &Scenario,
    spec: &NetworkSpec,
    steps: u64,
    dir: &Path,
    at: u64,
) -> Result<String> {
    let real = run(scenario.clone(), spec, steps)?;
    let lineage_of = |clone: CloneId| -> Vec<&CloneRecord> {
        let mut chain = Vec::new();
        let mut here = real.clones.get(clone.0 as usize);
        while let Some(record) = here {
            chain.push(record);
            here = record.parent.and_then(|p| real.clones.get(p.0 as usize));
        }
        chain
    };
    // Every mutation that was acting when the lesion appeared, in the lineage
    // of a clone that touched it.
    let candidates: BTreeSet<CloneId> = real
        .touched
        .iter()
        .flat_map(|&c| lineage_of(c))
        .filter(|r| r.mutation.is_some() && !r.suppressed && r.founded_at < at)
        .map(|r| r.id)
        .collect();
    // Cells carrying each candidate when the lesion appeared.
    let carrying = |candidate: CloneId| -> usize {
        real.sizes_at_lesion
            .iter()
            .enumerate()
            .filter(|&(c, _)| {
                lineage_of(CloneId(c as u32))
                    .iter()
                    .any(|r| r.id == candidate)
            })
            .map(|(_, n)| n)
            .sum()
    };
    let mut jobs: Vec<Vec<CloneId>> = candidates.iter().map(|&c| vec![c]).collect();
    if candidates.len() > 1 {
        jobs.push(candidates.iter().copied().collect());
    }
    let outcomes: Vec<LesionOutcome> = jobs
        .par_iter()
        .map(|undone| -> Result<LesionOutcome> {
            let mut reverted = scenario.clone();
            reverted.revert = Some(Reversion {
                clones: undone.clone(),
                at,
            });
            let other = run(reverted, spec, steps)?;
            let parted_at = real
                .measured
                .iter()
                .zip(&other.measured)
                .find(|(a, b)| a != b)
                .map(|(a, _)| a.mcs);
            let last = other
                .measured
                .last()
                .context("a run takes at least one step")?;
            Ok(LesionOutcome {
                undone: undone.clone(),
                identical: parted_at.is_none()
                    && real.measured.len() == other.measured.len()
                    && real.lattice == other.lattice,
                parted_at,
                lesion: last.lesion,
                lost_at: other.lesion_lost_at,
                end_mcs: last.mcs,
            })
        })
        .collect::<Result<_>>()?;

    let last = real
        .measured
        .last()
        .context("a run takes at least one step")?;
    let appeared = real
        .measured
        .iter()
        .find(|m| m.mcs == at)
        .map_or(0, |m| m.lesion);
    let fate = |lesion: usize, lost_at: Option<u64>, end: u64| match lost_at {
        Some(step) => format!("lost at step {step}"),
        None => format!("{lesion} cells at step {end}"),
    };
    let mut text = format!(
        "real history: the lesion appeared at step {at} over {appeared} cells and was {}; \
         {} clones touched it, carrying {} mutations made before it appeared\n",
        match real.lesion_lost_at {
            Some(step) => format!("lost at step {step}"),
            None => format!("still {} cells at step {}", last.lesion, last.mcs),
        },
        real.touched.len(),
        candidates.len()
    );
    let mut csv = String::from(
        "undone,founded_at,mutation,carrying,identical,parted_at,lesion_at_end,lesion_lost_at\n",
    );
    for o in &outcomes {
        let effect = if o.identical {
            "no effect: identical to the last bit".to_owned()
        } else {
            format!(
                "parted at step {}; the lesion {}",
                o.parted_at
                    .map_or_else(|| "?".to_owned(), |s| s.to_string()),
                fate(o.lesion, o.lost_at, o.end_mcs)
            )
        };
        let (label, founded, mutation, carried) = match o.undone.as_slice() {
            [one] => {
                let record = real.clones.get(one.0 as usize);
                (
                    format!("clone {}", one.0),
                    record.map_or(0, |r| r.founded_at).to_string(),
                    record
                        .and_then(|r| r.mutation.as_ref())
                        .map_or_else(String::new, ToString::to_string),
                    carrying(*one),
                )
            }
            all => (
                format!("all {}", all.len()),
                String::new(),
                String::new(),
                real.sizes_at_lesion
                    .iter()
                    .enumerate()
                    .filter(|&(c, _)| {
                        lineage_of(CloneId(c as u32))
                            .iter()
                            .any(|r| all.contains(&r.id))
                    })
                    .map(|(_, n)| n)
                    .sum(),
            ),
        };
        let what = if founded.is_empty() {
            format!("{label} at once ({carried} cells)")
        } else {
            format!("{label} (founded at step {founded}, {mutation}; {carried} cells)")
        };
        let _ = writeln!(text, "  undo {what}: {effect}");
        let _ = writeln!(
            csv,
            "{},{},{},{},{},{},{},{}",
            label,
            founded,
            mutation,
            carried,
            o.identical,
            o.parted_at.map_or_else(String::new, |s| s.to_string()),
            o.lesion,
            o.lost_at.map_or_else(String::new, |s| s.to_string())
        );
    }
    fs::create_dir_all(dir).with_context(|| format!("creating {}", dir.display()))?;
    fs::write(dir.join("attribution.csv"), csv).context("writing attribution.csv")?;
    Ok(text)
}
