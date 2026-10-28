//! One configuration under many seeds, or two configurations under the same
//! seeds, side by side.
//!
//! Two arms run under the same seed are the same history until the first
//! thing that differs between their configurations, say the step an exposure
//! goes off in one and stays on in the other. Comparing them seed by seed,
//! rather than mean against mean, takes that shared history out of the
//! comparison: what differs between the two runs of a pair is what differs
//! between the arms, and not the luck of the draw. That is the method of
//! common random numbers, and the exact reproducibility of a run is what
//! makes it available.

use std::fmt::Write as _;
use std::fs;
use std::ops::RangeInclusive;
use std::path::Path;

use anyhow::{Context, Result, bail};
use cellweave_core::NetworkSpec;
use cellweave_engine::simulation::Scenario;
use cellweave_engine::{DEPTH_SHARE, Simulation};
use rayon::prelude::*;

/// Fraction of the exposure's intensity a cell's response is measured at.
const RESPONSE_LEVEL: f64 = 0.5;

/// One arm of a campaign: a name and what it runs.
pub struct Arm {
    /// Shown in the table and the file.
    pub name: String,
    /// The configuration, under whatever seed the campaign sets.
    pub scenario: Scenario,
    /// The network its cells carry.
    pub spec: NetworkSpec,
    /// Steps per run.
    pub steps: u64,
}

/// Where one run ended.
#[derive(Debug, Clone, PartialEq)]
struct Outcome {
    seed: u64,
    arm: usize,
    end_mcs: u64,
    reached_edge: bool,
    living: usize,
    cycling: usize,
    quiescent: usize,
    clones: usize,
    /// Fewest living cells from the exposure's start on, when there is one.
    nadir: Option<usize>,
    /// Share of the living that is resistant at the end, when the network
    /// has an apoptosis role and the run has an exposure.
    resistant: Option<f64>,
    /// Mean response of the living to the stress at the end, relative to
    /// the unmutated network's: one unchanged, lower less sensitive.
    response: Option<f64>,
    /// The latest mutations of the clone holding the most living cells.
    leading: String,
    /// For each node, the share of the living whose lineage carries a change
    /// to it, as [`Simulation::coverage`] measures it.
    coverage: Vec<(String, f64)>,
    /// Living cells descended from the lesion at the end, when the
    /// configuration places one.
    lesion: Option<usize>,
    /// The first step at whose end no cell of the lesion was left, if it was
    /// lost. A lineage with no living cell cannot come back.
    lesion_lost_at: Option<u64>,
    /// Dead cells not yet resorbed, at the end.
    necrotic: usize,
    /// Pixels of vessel at the end, those of them carrying blood, and
    /// sprouts ever started, when the configuration has a vessel.
    vessel: Option<(u64, u64, usize)>,
    /// How far living tissue reaches from its nearest vessel at the end, as
    /// [`Simulation::depth_from_vessels`] reads it.
    depth: Option<f64>,
}

/// Mutations shown for the leading clone, the latest ones.
const SHOWN_MUTATIONS: usize = 4;

/// A range of seeds written `a-b`, or a single seed.
pub fn parse_seeds(text: &str) -> Result<RangeInclusive<u64>> {
    let (a, b) = text.split_once('-').unwrap_or((text, text));
    let (a, b): (u64, u64) = (
        a.trim()
            .parse()
            .context("--seeds takes a seed or a range like 1-20")?,
        b.trim()
            .parse()
            .context("--seeds takes a seed or a range like 1-20")?,
    );
    if a > b {
        bail!("--seeds {a}-{b} is empty; write the smaller seed first");
    }
    Ok(a..=b)
}

fn run(arm: &Arm, index: usize, seed: u64) -> Result<Outcome> {
    let scenario = Scenario {
        seed,
        ..arm.scenario.clone()
    };
    let exposure = scenario.exposure;
    let lesion_at = scenario.lesion.as_ref().map(|l| l.at);
    let (mut sim, _) = Simulation::new(scenario, &arm.spec)?.start()?;
    let mut lesion_lost_at = None;
    let mut nadir: Option<usize> = None;
    let mut reached_edge = false;
    for _ in 0..arm.steps {
        let step = sim.step()?;
        if exposure.is_some_and(|e| step.mcs >= e.start) {
            let living = sim.living();
            nadir = Some(nadir.map_or(living, |n| n.min(living)));
        }
        if lesion_lost_at.is_none()
            && lesion_at.is_some_and(|at| step.mcs >= at)
            && sim.lesion_cells() == 0
        {
            lesion_lost_at = Some(step.mcs);
        }
        if step.reached_edge {
            reached_edge = true;
            break;
        }
    }
    let m = sim.measure();
    // Asked at half the source's level, which is about what the cells nearest
    // the source read once the stress has been taken up on its way in; at the
    // full level a response saturates, and no single mutation can halve it.
    let level = exposure.map(|e| RESPONSE_LEVEL * e.intensity);
    let resistant = level.and_then(|l| sim.resistant_fraction(l));
    let response = level
        .and_then(|l| sim.relative_responses(l))
        .filter(|r| !r.is_empty())
        .map(|r| r.iter().sum::<f64>() / r.len() as f64);
    let sizes = sim.clone_sizes();
    let leading = (0..sizes.len())
        .filter(|&c| sizes[c] > 0)
        .max_by_key(|&c| (sizes[c], std::cmp::Reverse(c)))
        .map_or_else(String::new, |c| {
            let mut chain = Vec::new();
            let mut at = sim.clones().get(c);
            while let Some(record) = at {
                if let Some(m) = &record.mutation
                    && !record.suppressed
                {
                    chain.push(m.to_string());
                }
                at = record.parent.and_then(|p| sim.clones().get(p.0 as usize));
            }
            chain.truncate(SHOWN_MUTATIONS);
            chain.reverse();
            chain.join(" / ")
        });
    Ok(Outcome {
        seed,
        arm: index,
        end_mcs: m.mcs,
        reached_edge,
        living: m.living(),
        cycling: m.cycling,
        quiescent: m.quiescent,
        clones: m.clones,
        nadir,
        resistant,
        response,
        leading,
        coverage: sim.coverage(),
        lesion: lesion_at.map(|_| m.lesion),
        lesion_lost_at,
        necrotic: m.necrotic,
        vessel: (m.vessel > 0).then_some((m.vessel, m.perfused, sim.sprouted())),
        depth: sim.depth_from_vessels(DEPTH_SHARE),
    })
}

/// Run every arm under every seed, all at once, write `campaign.csv` into
/// `dir`, and return a summary: each arm on its own, and with two arms, how
/// the second compares with the first seed by seed.
pub fn campaign(arms: &[Arm], seeds: RangeInclusive<u64>, dir: &Path) -> Result<String> {
    let jobs: Vec<(usize, u64)> = seeds
        .clone()
        .flat_map(|s| (0..arms.len()).map(move |a| (a, s)))
        .collect();
    let mut outcomes = jobs
        .par_iter()
        .map(|&(a, s)| run(&arms[a], a, s))
        .collect::<Result<Vec<_>>>()?;
    outcomes.sort_by_key(|o| (o.seed, o.arm));

    let mut csv = String::from(
        "seed,arm,end_mcs,reached_edge,living,cycling,quiescent,clones,nadir,resistant,response,leading,lesion,lesion_lost_at,necrotic,vessel,perfused,sprouted,depth,coverage\n",
    );
    for o in &outcomes {
        let _ = writeln!(
            csv,
            "{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{},{}",
            o.seed,
            arms[o.arm].name,
            o.end_mcs,
            o.reached_edge,
            o.living,
            o.cycling,
            o.quiescent,
            o.clones,
            o.nadir.map_or_else(String::new, |n| n.to_string()),
            o.resistant.map_or_else(String::new, |r| format!("{r:.4}")),
            o.response.map_or_else(String::new, |r| format!("{r:.4}")),
            o.leading,
            o.lesion.map_or_else(String::new, |n| n.to_string()),
            o.lesion_lost_at.map_or_else(String::new, |s| s.to_string()),
            o.necrotic,
            o.vessel.map_or_else(String::new, |(v, _, _)| v.to_string()),
            o.vessel.map_or_else(String::new, |(_, p, _)| p.to_string()),
            o.vessel.map_or_else(String::new, |(_, _, s)| s.to_string()),
            o.depth.map_or_else(String::new, |d| format!("{d:.1}")),
            o.coverage
                .iter()
                .filter(|(_, share)| *share > 0.0)
                .map(|(node, share)| format!("{node}={share:.4}"))
                .collect::<Vec<_>>()
                .join(" ")
        );
    }
    fs::create_dir_all(dir).with_context(|| format!("creating {}", dir.display()))?;
    fs::write(dir.join("campaign.csv"), csv).context("writing campaign.csv")?;

    let n = seeds.clone().count();
    let mut text = format!("{n} seeds, {}..={}\n", seeds.start(), seeds.end());
    for (a, arm) in arms.iter().enumerate() {
        let mine: Vec<&Outcome> = outcomes.iter().filter(|o| o.arm == a).collect();
        let living: Vec<f64> = mine.iter().map(|o| o.living as f64).collect();
        let edge = mine.iter().filter(|o| o.reached_edge).count();
        let resistant: Vec<f64> = mine.iter().filter_map(|o| o.resistant).collect();
        let response: Vec<f64> = mine.iter().filter_map(|o| o.response).collect();
        let _ = write!(
            text,
            "  {}: living at the end {}, touched the edge in {edge} of {n}",
            arm.name,
            spread(&living)
        );
        if !resistant.is_empty() {
            let _ = write!(
                text,
                ", resistant share {}, response {}",
                spread(&resistant),
                spread(&response)
            );
        }
        text.push('\n');
        let vessel: Vec<f64> = mine
            .iter()
            .filter_map(|o| o.vessel)
            .map(|(v, _, _)| v as f64)
            .collect();
        if !vessel.is_empty() {
            let perfused: Vec<f64> = mine
                .iter()
                .filter_map(|o| o.vessel)
                .map(|(v, p, _)| p as f64 / v as f64)
                .collect();
            let sprouted: Vec<f64> = mine
                .iter()
                .filter_map(|o| o.vessel)
                .map(|(_, _, s)| s as f64)
                .collect();
            let necrotic: Vec<f64> = mine.iter().map(|o| o.necrotic as f64).collect();
            let depth: Vec<f64> = mine.iter().filter_map(|o| o.depth).collect();
            let _ = writeln!(
                text,
                "    necrotic {}, vessel pixels {}, carrying blood {}, sprouts {}, living tissue within {} px of a vessel",
                spread(&necrotic),
                spread(&vessel),
                spread(&perfused),
                spread(&sprouted),
                spread(&depth)
            );
        }
        if let Some(at) = arm.scenario.lesion.as_ref().map(|l| l.at) {
            let kept = mine
                .iter()
                .filter(|o| o.lesion.is_some_and(|n| n > 0))
                .count();
            let lost: Vec<f64> = mine
                .iter()
                .filter_map(|o| o.lesion_lost_at)
                .map(|s| (s - at) as f64)
                .collect();
            let _ = write!(
                text,
                "    lesion: still there at the end in {kept} of {n}, lost in {}",
                lost.len()
            );
            if !lost.is_empty() {
                let _ = write!(text, ", steps after it appeared: {}", spread(&lost));
            }
            text.push('\n');
        }
        // Every node any run of this arm carries a change to, with the
        // share of the living that carries it, seed by seed.
        let nodes: Vec<&str> = mine
            .iter()
            .flat_map(|o| o.coverage.iter())
            .filter(|(_, share)| *share > 0.0)
            .map(|(node, _)| node.as_str())
            .fold(Vec::new(), |mut seen, node| {
                if !seen.contains(&node) {
                    seen.push(node);
                }
                seen
            });
        for node in nodes {
            let shares: Vec<f64> = mine
                .iter()
                .map(|o| {
                    o.coverage
                        .iter()
                        .find(|(n, _)| n == node)
                        .map_or(0.0, |(_, share)| *share)
                })
                .collect();
            let _ = writeln!(text, "    carrying a change to {node}: {}", spread(&shares));
        }
    }
    if arms.len() == 2 {
        let pairs: Vec<(&Outcome, &Outcome)> = outcomes
            .chunks(2)
            .filter_map(|p| match p {
                [a, b] if a.seed == b.seed => Some((a, b)),
                _ => None,
            })
            .collect();
        let (first, second) = (&arms[0].name, &arms[1].name);
        let _ = writeln!(text, "  seed by seed, {second} against {first}:");
        let _ = writeln!(
            text,
            "    more living at the end: {}",
            tally(&pairs, |o| Some(o.living as f64))
        );
        let _ = writeln!(
            text,
            "    larger resistant share: {}",
            tally(&pairs, |o| o.resistant)
        );
        let _ = writeln!(
            text,
            "    weaker response, more resistant: {}",
            tally(&pairs, |o| o.response.map(|r| -r))
        );
        let _ = writeln!(
            text,
            "    more lesion cells at the end: {}",
            tally(&pairs, |o| o.lesion.map(|n| n as f64))
        );
        if pairs
            .iter()
            .all(|(a, b)| a.depth.is_some() && b.depth.is_some())
        {
            let _ = writeln!(
                text,
                "    more necrotic at the end: {}",
                tally(&pairs, |o| Some(o.necrotic as f64))
            );
            let _ = writeln!(
                text,
                "    living tissue reaching further from a vessel: {}",
                tally(&pairs, |o| o.depth)
            );
        }
    }
    Ok(text)
}

/// Median and range of some numbers, as text.
fn spread(values: &[f64]) -> String {
    if values.is_empty() {
        return "none".to_owned();
    }
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    let median = if sorted.len() % 2 == 1 {
        sorted[sorted.len() / 2]
    } else {
        (sorted[sorted.len() / 2 - 1] + sorted[sorted.len() / 2]) / 2.0
    };
    format!(
        "median {:.2} (from {:.2} to {:.2})",
        median,
        sorted[0],
        sorted[sorted.len() - 1]
    )
}

/// In how many pairs the second arm's value is above, below, or equal to
/// the first's, counting only pairs where both have one.
fn tally(pairs: &[(&Outcome, &Outcome)], pick: impl Fn(&Outcome) -> Option<f64>) -> String {
    let (mut above, mut below, mut equal) = (0, 0, 0);
    for (a, b) in pairs {
        if let (Some(x), Some(y)) = (pick(a), pick(b)) {
            match y.total_cmp(&x) {
                std::cmp::Ordering::Greater => above += 1,
                std::cmp::Ordering::Less => below += 1,
                std::cmp::Ordering::Equal => equal += 1,
            }
        }
    }
    format!("in {above} seeds, fewer in {below}, equal in {equal}")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn seeds_are_a_seed_or_a_range() {
        assert_eq!(parse_seeds("7").unwrap(), 7..=7);
        assert_eq!(parse_seeds("1-20").unwrap(), 1..=20);
        assert!(parse_seeds("20-1").is_err());
        assert!(parse_seeds("one").is_err());
    }

    #[test]
    fn a_spread_is_a_median_and_a_range() {
        assert_eq!(spread(&[3.0, 1.0, 2.0]), "median 2.00 (from 1.00 to 3.00)");
        assert_eq!(spread(&[4.0, 1.0]), "median 2.50 (from 1.00 to 4.00)");
        assert_eq!(spread(&[]), "none");
    }
}
