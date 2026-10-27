//! A lethal stress and apoptosis, end to end: the shipped network and the
//! configuration it comes with, run on a small bath.

#![allow(clippy::unwrap_used)]

use cellweave_core::{
    CellKind, DiffusionField, Exposure, NetworkSpec, OxygenSource, Schedule, SignalingNetwork,
};
use cellweave_engine::simulation::{DEFAULT_SEED, Running, Scenario};
use cellweave_engine::{Fate, Simulation};
use cellweave_io::SimConfig;

const CONFIG: &str = concat!(
    env!("CARGO_MANIFEST_DIR"),
    "/../../examples/cord-stress.toml"
);

fn spec() -> NetworkSpec {
    let config = SimConfig::from_file(CONFIG).unwrap();
    config.load_network(Some(CONFIG)).unwrap()
}

fn exposed(exposure: Option<Exposure>) -> Scenario {
    Scenario {
        width: 80,
        height: 80,
        temperature: 10.0,
        oxygen_source: OxygenSource::Bath,
        oxygen_uptake: 0.0003,
        seed: DEFAULT_SEED,
        mutation_rate: 0.0,
        suppress: Vec::new(),
        rewiring: 0.0,
        seeding: Default::default(),
        exposure,
    }
}

fn continuous(start: u64) -> Option<Exposure> {
    Some(Exposure {
        intensity: 1.0,
        uptake: 0.001,
        start,
        schedule: Schedule::Continuous,
    })
}

fn run(scenario: Scenario, steps: u32) -> Simulation<Running> {
    let (mut sim, _) = Simulation::new(scenario, &spec()).unwrap().start().unwrap();
    for _ in 0..steps {
        sim.step().unwrap();
    }
    sim
}

#[test]
fn a_schedule_is_on_when_it_says() {
    let periodic = Exposure {
        intensity: 0.8,
        uptake: 0.001,
        start: 10,
        schedule: Schedule::Periodic { on: 3, off: 2 },
    };
    let levels: Vec<f64> = (8..20).map(|mcs| periodic.level_at(mcs)).collect();
    #[rustfmt::skip]
    let expected = [0.0, 0.0, 0.8, 0.8, 0.8, 0.0, 0.0, 0.8, 0.8, 0.8, 0.0, 0.0];
    assert_eq!(levels, expected);
    assert_eq!(continuous(5).unwrap().level_at(4), 0.0);
    assert_eq!(continuous(5).unwrap().level_at(500), 1.0);
}

#[test]
fn without_exposure_there_is_no_field_and_no_apoptosis() {
    let sim = run(exposed(None), 120);
    assert!(sim.stress().is_none());
    assert_eq!(sim.measure().apoptotic, 0);
    assert_eq!(sim.measure().exposure, 0.0);
}

#[test]
fn the_stress_is_taken_up_on_its_way_in_and_kills_by_apoptosis() {
    // The tissue grows for a hundred steps in peace, then is exposed. The
    // field is looked at two steps in, while every cell still takes it up.
    let mut sim = run(exposed(continuous(100)), 102);
    let field = sim.stress().unwrap();
    let (w, h) = field.dims();
    // Held full at the edge of the box, where the bath is, and lower inside
    // the tissue, which has taken some of it up on the way.
    assert!((field.concentration_at(0, 0) - 1.0).abs() < 1e-9);
    let lowest = (0..h)
        .flat_map(|y| (0..w).map(move |x| (x, y)))
        .map(|(x, y)| field.concentration_at(x, y))
        .fold(f64::INFINITY, f64::min);
    assert!(lowest < 0.99, "the tissue took some up: lowest {lowest}");

    for _ in 0..58 {
        sim.step().unwrap();
    }

    let lattice = sim.lattice();
    let apoptotic = lattice
        .cell_records()
        .filter(|(_, r)| r.kind == CellKind::Apoptotic)
        .count();
    assert!(
        apoptotic > 0 || sim.gone() > 0,
        "some cell started apoptosis under the stress"
    );
    assert_eq!(sim.measure().exposure, 1.0);
}

#[test]
fn an_apoptotic_cell_breaks_up_and_is_gone() {
    let sim = run(exposed(continuous(100)), 260);
    // Nothing dies of anything else in a bath this small this early, so what
    // is gone was cleared after apoptosis.
    assert_eq!(sim.necrotic(), 0);
    assert!(sim.gone() > 0, "apoptotic cells are cleared");
    let breaking_up = (1..=sim.lattice().cell_count() as u32)
        .map(cellweave_core::CellId)
        .filter(|&id| sim.fate(id) == Some(Fate::Apoptotic))
        .count();
    assert_eq!(sim.measure().apoptotic, breaking_up);
}

#[test]
fn an_exposure_that_cannot_settle_is_refused() {
    let bad = |e: Exposure| Simulation::new(exposed(Some(e)), &spec()).is_err();
    let base = continuous(0).unwrap();
    assert!(
        bad(Exposure {
            uptake: 0.0,
            ..base
        }),
        "nothing takes it up"
    );
    assert!(
        bad(Exposure {
            intensity: -1.0,
            ..base
        }),
        "a negative level"
    );
    assert!(
        bad(Exposure {
            intensity: f64::NAN,
            ..base
        }),
        "not a number"
    );
    assert!(bad(Exposure {
        schedule: Schedule::Periodic { on: 0, off: 5 },
        ..base
    }));
    let field_off = Scenario {
        oxygen_uptake: 0.0,
        ..exposed(continuous(0))
    };
    assert!(
        Simulation::new(field_off, &spec()).is_err(),
        "no source to come in from"
    );
}

#[test]
fn a_response_is_read_on_a_copy_and_needs_the_role() {
    let net = cellweave_engine::Network::from_spec(&spec(), 200).unwrap();
    let inputs = cellweave_core::CellInputs {
        o2: cellweave_core::Concentration(1.0),
        egf: cellweave_core::Concentration(0.2),
        vegf: cellweave_core::Concentration(0.0),
        crowding: cellweave_core::Concentration(0.0),
        stress: cellweave_core::Concentration(0.5),
    };
    let before = net.state().values.clone();
    let chance = net.apoptosis_under(&inputs, 60).unwrap();
    assert!(
        chance > 0.05,
        "the unmutated network answers half the source with {chance}"
    );
    assert_eq!(
        net.state().values,
        before,
        "asking does not move the network"
    );

    let plain: NetworkSpec =
        toml::from_str(include_str!("../../../examples/networks/ras-erk.toml")).unwrap();
    let plain = cellweave_engine::Network::from_spec(&plain, 200).unwrap();
    assert_eq!(plain.apoptosis_under(&inputs, 60), None);
}

#[test]
fn turning_p53_down_makes_a_cell_resistant() {
    let net = cellweave_engine::Network::from_spec(&spec(), 200).unwrap();
    let inputs = cellweave_core::CellInputs {
        o2: cellweave_core::Concentration(1.0),
        egf: cellweave_core::Concentration(0.2),
        vegf: cellweave_core::Concentration(0.0),
        crowding: cellweave_core::Concentration(0.0),
        stress: cellweave_core::Concentration(0.5),
    };
    // Sixteen terms, the twelve of the RAS/ERK network and then p53's two and
    // caspase's two: a draw of 12.5 / 16 lands on p53's activation, and a draw
    // of zero halves it.
    let mut mutant = net.clone();
    let m = mutant.mutate(12.5 / 16.0, 0.0).unwrap();
    assert_eq!((m.node.as_str(), m.kind), ("p53", "activated_by"));
    let ratio =
        mutant.apoptosis_under(&inputs, 60).unwrap() / net.apoptosis_under(&inputs, 60).unwrap();
    assert!(
        ratio < cellweave_engine::evolution::RESISTANT_BELOW,
        "half the p53 makes a cell resistant: {ratio}"
    );
}

#[test]
fn an_unmutated_population_has_no_resistant_share() {
    let sim = run(exposed(continuous(100)), 101);
    assert_eq!(sim.resistant_fraction(0.5), Some(0.0));
    let plain: NetworkSpec =
        toml::from_str(include_str!("../../../examples/networks/ras-erk.toml")).unwrap();
    let (without_role, _) = Simulation::new(exposed(None), &plain)
        .unwrap()
        .start()
        .unwrap();
    assert_eq!(without_role.resistant_fraction(0.5), None);
}
