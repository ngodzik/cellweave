//! A sheet of normal epithelium, end to end: read from the shipped files, run,
//! and checked for what an epithelium has to do.

#![allow(clippy::unwrap_used)]

use cellweave_core::{CellId, CloneId, NetworkSpec, OxygenSource, Seeding, Tracing};
use cellweave_engine::simulation::{DEFAULT_SEED, Running, Scenario};
use cellweave_engine::{Fate, Simulation};
use cellweave_io::SimConfig;

const CONFIG: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../examples/epithelium.toml"
);

fn spec() -> NetworkSpec {
    let config = SimConfig::from_file(CONFIG).unwrap();
    config.load_network(Some(CONFIG)).unwrap()
}

fn sheet(size: u32) -> Scenario {
    Scenario {
        width: size,
        height: size,
        temperature: 10.0,
        oxygen_source: OxygenSource::Bath,
        oxygen_uptake: 0.0,
        seed: DEFAULT_SEED,
        mutation_rate: 0.0,
        suppress: Vec::new(),
        suppress_before: None,
        rewiring: 0.0,
        seeding: Seeding::Sheet { side: 10 },
        exposure: None,
        lesion: None,
        revert: None,
        angiogenesis: None,
        immunity: None,
        tracing: None,
    }
}

fn run(size: u32, steps: u32) -> Simulation<Running> {
    let (mut sim, sweeps) = Simulation::new(sheet(size), &spec())
        .unwrap()
        .start()
        .unwrap();
    assert_eq!(sweeps, 0, "nothing consumes, so the field is not solved");
    for _ in 0..steps {
        let step = sim.step().unwrap();
        assert!(
            !step.reached_edge,
            "a sheet starts at the edge and is not stopped by it"
        );
    }
    sim
}

#[test]
fn a_sheet_fills_the_box_with_one_clone_per_cell() {
    let sim = Simulation::new(sheet(40), &spec()).unwrap();
    assert_eq!(sim.lattice().cell_count(), 16);
    for y in 0..40 {
        for x in 0..40 {
            assert_ne!(sim.lattice().occupant(x, y), 0, "no medium at ({x}, {y})");
        }
    }
    for i in 1..=16 {
        assert_eq!(sim.clone_of(CellId(i)), Some(CloneId(i - 1)));
    }
    assert_eq!(sim.clones().len(), 16);
}

#[test]
fn a_sheet_that_consumes_oxygen_is_refused() {
    let scenario = Scenario {
        oxygen_uptake: 0.0003,
        ..sheet(40)
    };
    assert!(Simulation::new(scenario, &spec()).is_err());
}

#[test]
fn the_layer_turns_over_and_holds_its_size() {
    let sim = run(60, 600);
    let start = 36;
    let living = sim.living();
    assert!(
        (start * 6 / 10..=start * 14 / 10).contains(&living),
        "{living} cells in the layer from {start}"
    );
    assert!(
        sim.born() > 3 * start as u64,
        "the layer turned over several times: {} born",
        sim.born()
    );
    assert!(sim.gone() > 0, "cells left the layer");
}

#[test]
fn with_no_advantage_clones_are_lost_by_chance_alone() {
    let sim = run(60, 600);
    let surviving = sim.clone_sizes().iter().filter(|&&n| n > 0).count();
    assert!(
        surviving < 36 / 2,
        "{surviving} of 36 clones still in the layer"
    );
    assert!(surviving >= 1);
}

#[test]
fn a_cell_that_commits_leaves_the_layer() {
    // Cells commit by chance, mostly just after a division: step until one is
    // on its way out rather than look at one step chosen in advance.
    let (mut sim, _) = Simulation::new(sheet(60), &spec())
        .unwrap()
        .start()
        .unwrap();
    let mut leaving: Vec<CellId> = Vec::new();
    for _ in 0..1000 {
        sim.step().unwrap();
        leaving = (1..=sim.lattice().cell_count() as u32)
            .map(CellId)
            .filter(|&id| sim.fate(id) == Some(Fate::Differentiating))
            .collect();
        if !leaving.is_empty() {
            break;
        }
    }
    assert!(
        !leaving.is_empty(),
        "some cell went on its way out within 1000 steps"
    );
    let measured = sim.measure();
    assert_eq!(measured.differentiating, leaving.len());
    assert_eq!(measured.quiescent, 0, "full oxygen arrests nothing");
}

#[test]
fn a_sheet_that_does_not_tile_the_box_is_refused() {
    // 64 is not a multiple of 10: the last row and column would be clipped
    // cells, crowded from the first step.
    assert!(Simulation::new(sheet(64), &spec()).is_err());
}

#[test]
fn a_sheet_around_a_vessel_is_refused() {
    let scenario = Scenario {
        oxygen_source: OxygenSource::Vessel { radius: 6 },
        ..sheet(40)
    };
    assert!(Simulation::new(scenario, &spec()).is_err());
}

#[test]
fn sisters_do_opposite_things_more_often_than_chance() {
    // Independent sisters, whatever their odds of dividing, give a quarter
    // or more: both divide or both leave at least half the time. Nothing in
    // the model couples two sisters but the room they share.
    let sim = run(200, 600);
    let fates = sim.sister_fates();
    assert!(fates.pairs() > 1000, "{fates:?}");
    let r = fates.r().unwrap();
    assert!(r < 0.24, "r = {r}");
}

#[test]
fn marking_gives_every_cell_a_clone_of_its_own_and_changes_nothing_else() {
    let traced = Scenario {
        tracing: Some(Tracing { at: 100, every: 10 }),
        ..sheet(60)
    };
    let (mut marked, _) = Simulation::new(traced, &spec()).unwrap().start().unwrap();
    let (mut plain, _) = Simulation::new(sheet(60), &spec())
        .unwrap()
        .start()
        .unwrap();
    let in_layer = |sim: &Simulation<Running>| {
        let m = sim.measure();
        m.cycling + m.quiescent + m.differentiating
    };
    for _ in 0..100 {
        marked.step().unwrap();
        plain.step().unwrap();
    }
    assert!(
        marked.traced_clones().is_empty(),
        "nothing is marked before the step"
    );

    marked.step().unwrap();
    plain.step().unwrap();
    // Marked at the start of the step, a cell may divide before its end.
    let at_marking = marked.traced_clones();
    assert!(at_marking.iter().all(|&(_, cells)| cells <= 2));
    assert_eq!(
        at_marking.iter().map(|&(_, cells)| cells).sum::<usize>(),
        in_layer(&marked),
        "every cell in the layer carries a mark"
    );

    for _ in 0..300 {
        marked.step().unwrap();
        plain.step().unwrap();
    }
    assert_eq!(
        marked.measure(),
        plain.measure(),
        "marking draws nothing, so the run is the unmarked run"
    );
    let later = marked.traced_clones();
    assert_eq!(
        later.iter().map(|&(_, cells)| cells).sum::<usize>(),
        in_layer(&marked),
        "every cell in the layer descends from a cell marked then"
    );
    assert!(
        later.len() < at_marking.len(),
        "marked clones are lost by drift"
    );
}

#[test]
fn tracing_every_zero_steps_is_refused() {
    let traced = Scenario {
        tracing: Some(Tracing { at: 10, every: 0 }),
        ..sheet(60)
    };
    assert!(Simulation::new(traced, &spec()).is_err());
}
