//! An early tumour in a sheet, end to end: read from the shipped files, placed,
//! followed, and checked for what placing it must and must not do.

#![allow(clippy::unwrap_used)]

use cellweave_core::{CellId, OxygenSource, Seeding};
use cellweave_engine::Simulation;
use cellweave_engine::simulation::{Lesion, Scenario};
use cellweave_io::SimConfig;

const CONFIG: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../examples/epithelium-lesion.toml"
);

/// The shipped lesion scenario on a smaller box, the lesion at `at`.
fn scenario(at: u64) -> (Scenario, cellweave_core::NetworkSpec) {
    let config = SimConfig::from_file(CONFIG).unwrap();
    let mut scenario = config.scenario(Some(CONFIG)).unwrap();
    scenario.width = 100;
    scenario.height = 100;
    if let Some(lesion) = scenario.lesion.as_mut() {
        lesion.at = at;
    }
    (scenario, config.load_network(Some(CONFIG)).unwrap())
}

#[test]
fn a_lesion_turns_the_cells_at_the_centre_into_one_clone_at_its_step() {
    let (scenario, spec) = scenario(20);
    let (mut sim, _) = Simulation::new(scenario, &spec).unwrap().start().unwrap();
    for _ in 0..20 {
        sim.step().unwrap();
    }
    assert_eq!(sim.lesion_cells(), 0, "not before its step");
    assert!(sim.clones().iter().all(|c| !c.lesion));
    sim.step().unwrap();
    let founded: Vec<_> = sim.clones().iter().filter(|c| c.lesion).collect();
    assert_eq!(founded.len(), 1, "one clone for the whole lesion");
    assert_eq!(founded[0].founded_at, 20);
    assert_eq!(founded[0].parent, None);
    let cells = sim.lesion_cells();
    assert!((3..=12).contains(&cells), "{cells} cells within 15 pixels");
    // The cell holding the centre of the box is one of them.
    let centre = CellId(sim.lattice().occupant(50, 50));
    assert!(sim.in_lesion(centre));
    assert!(!sim.in_lesion(CellId(sim.lattice().occupant(2, 2))));
}

#[test]
fn a_lesion_that_never_comes_leaves_the_run_as_it_was() {
    let run = |lesion: bool| {
        let (mut scenario, spec) = scenario(10_000);
        if !lesion {
            scenario.lesion = None;
        }
        let (mut sim, _) = Simulation::new(scenario, &spec).unwrap().start().unwrap();
        let mut measured = Vec::new();
        for _ in 0..200 {
            sim.step().unwrap();
            measured.push(sim.measure());
        }
        (measured, sim.clones().to_vec())
    };
    assert_eq!(run(true), run(false), "placing a lesion draws nothing");
}

#[test]
fn a_lesion_needs_a_sheet_and_some_room() {
    let (scenario, spec) = scenario(20);
    let single = Scenario {
        seeding: Seeding::Single,
        oxygen_source: OxygenSource::Bath,
        oxygen_uptake: 0.0003,
        ..scenario.clone()
    };
    assert!(Simulation::new(single, &spec).is_err());
    let empty = Scenario {
        lesion: scenario.lesion.clone().map(|l| Lesion { radius: 0.0, ..l }),
        ..scenario
    };
    assert!(Simulation::new(empty, &spec).is_err());
}
