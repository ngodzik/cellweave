//! An early tumour in a sheet, end to end: read from the shipped files, placed,
//! followed, and checked for what placing it must and must not do; and a
//! patch placed by hand in a tumour, off centre.

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
fn a_lesion_needs_some_room_inside_the_box() {
    let (scenario, spec) = scenario(20);
    let empty = Scenario {
        lesion: scenario.lesion.clone().map(|l| Lesion { radius: 0.0, ..l }),
        ..scenario.clone()
    };
    assert!(Simulation::new(empty, &spec).is_err());
    let outside = Scenario {
        lesion: scenario.lesion.clone().map(|l| Lesion {
            centre: Some((150.0, 50.0)),
            ..l
        }),
        ..scenario
    };
    assert!(Simulation::new(outside, &spec).is_err());
}

#[test]
fn a_lesion_placed_off_centre_in_a_tumour_takes_the_cells_there() {
    let config = SimConfig::default();
    let spec = config.load_network(None).unwrap();
    // One cell seeded against a vessel of radius 6, eight pixels out from
    // its wall: the centre of the box is vessel, and the cell is to its right.
    let first = (100.0 + 6.0 + 8.0, 100.0);
    let tumour = Scenario {
        seeding: Seeding::Single,
        oxygen_source: OxygenSource::Vessel { radius: 6 },
        oxygen_uptake: 0.00003,
        lesion: Some(Lesion {
            spec: config.load_network(None).unwrap(),
            at: 1,
            radius: 4.0,
            centre: Some(first),
        }),
        ..config.scenario(None).unwrap()
    };
    let (mut sim, _) = Simulation::new(tumour, &spec).unwrap().start().unwrap();
    sim.step().unwrap();
    assert_eq!(sim.lesion_cells(), 0, "not before its step");
    sim.step().unwrap();
    assert_eq!(sim.lesion_cells(), 1, "the one cell there");
    assert!(sim.in_lesion(CellId(1)));
}
