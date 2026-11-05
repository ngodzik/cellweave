//! The fuller cell of chapter 28 against the networks it fills in: the same
//! decisions, at the same places, with more proteins between.

#![allow(clippy::unwrap_used)]

use cellweave_core::{CellInputs, Concentration, NetworkSpec, SignalingNetwork, TimeStep};
use cellweave_engine::Network;
use cellweave_engine::simulation::QUIESCENCE_THRESHOLD;
use cellweave_io::SimConfig;

fn spec(config: &str) -> NetworkSpec {
    let path = format!("{}/../../examples/{config}", env!("CARGO_MANIFEST_DIR"));
    let config = SimConfig::from_file(&path).unwrap();
    config.load_network(Some(&path)).unwrap()
}

fn fuller() -> NetworkSpec {
    spec("spheroid-fuller.toml")
}

fn simpler() -> NetworkSpec {
    spec("spheroid.toml")
}

fn stress_response() -> NetworkSpec {
    spec("cord-stress.toml")
}

fn inputs(o2: f64, stress: f64) -> CellInputs {
    CellInputs {
        o2: Concentration(o2),
        egf: Concentration(0.2),
        vegf: Concentration(0.0),
        crowding: Concentration(0.0),
        stress: Concentration(stress),
        waste: Concentration(0.0),
        necrotic_factor: Concentration(0.0),
    }
}

/// The network of `spec`, settled under `o2` and `stress`.
fn settled(spec: &NetworkSpec, o2: f64, stress: f64) -> Network {
    let mut net = Network::from_spec(spec, 100).unwrap();
    for _ in 0..600 {
        net.step(TimeStep(1.0), &inputs(o2, stress));
    }
    net
}

fn level(net: &Network, name: &str) -> f64 {
    let i = net.node_names().iter().position(|n| n == name).unwrap();
    net.state().values[i].0
}

/// The oxygen below which a settled cell of `spec` has stopped dividing.
fn arrest_oxygen(spec: &NetworkSpec) -> f64 {
    let (mut lo, mut hi) = (0.0, 1.0);
    for _ in 0..30 {
        let mid = 0.5 * (lo + hi);
        if settled(spec, mid, 0.0).hypoxia_response() >= QUIESCENCE_THRESHOLD {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

#[test]
fn the_fuller_cell_asks_for_the_size_the_simpler_one_does() {
    let (a, b) = (settled(&fuller(), 1.0, 0.0), settled(&simpler(), 1.0, 0.0));
    assert!(
        (level(&a, "erk") - level(&b, "erk")).abs() < 0.002,
        "ERK {} against {}",
        level(&a, "erk"),
        level(&b, "erk")
    );
    let (ta, tb) = (a.mechanics().target_volume.0, b.mechanics().target_volume.0);
    assert!(ta.abs_diff(tb) <= 1, "target {ta} against {tb}");
}

#[test]
fn the_fuller_cell_stops_dividing_at_the_oxygen_the_simpler_one_does() {
    let (a, b) = (arrest_oxygen(&fuller()), arrest_oxygen(&simpler()));
    assert!((a - b).abs() < 0.005, "arrest below {a:.4} against {b:.4}");
    // And by RB, through p27: well fed, RB is off and p27 absent.
    let fed = settled(&fuller(), 1.0, 0.0);
    assert!(level(&fed, "rb") < 0.3, "RB {}", level(&fed, "rb"));
    assert!(level(&fed, "p27") < 0.01);
    let starved = settled(&fuller(), 0.1, 0.0);
    assert!(
        level(&starved, "p27") > 0.5,
        "p27 {}",
        level(&starved, "p27")
    );
    assert!(level(&starved, "cyclin_e") < level(&fed, "cyclin_e") / 2.0);
}

#[test]
fn the_fuller_cell_dies_of_anoxia_where_the_simpler_one_does() {
    for o2 in [0.02, 0.05, 0.1, 0.2, 0.6] {
        let (a, b) = (settled(&fuller(), o2, 0.0), settled(&simpler(), o2, 0.0));
        assert!(
            (a.survival() - b.survival()).abs() < 0.002,
            "at {o2}: BCL2 {} against {}",
            a.survival(),
            b.survival()
        );
    }
}

#[test]
fn the_fuller_cell_answers_a_stress_as_the_stress_response_did() {
    for stress in [0.0, 0.2, 0.5, 1.0] {
        let (a, b) = (
            settled(&fuller(), 1.0, stress),
            settled(&stress_response(), 1.0, stress),
        );
        let (pa, pb) = (a.apoptosis().unwrap(), b.apoptosis().unwrap());
        assert!(
            (pa - pb).abs() < 0.015,
            "at {stress}: caspase {pa:.4} against {pb:.4}"
        );
    }
    // Without a stress, nothing on the way to apoptosis is on at all.
    let calm = settled(&fuller(), 1.0, 0.0);
    for name in ["p53", "mdm2", "p21", "bax", "caspase"] {
        assert_eq!(level(&calm, name), 0.0, "{name}");
    }
    // MDM2 holds p53 below what the stress alone would raise it to.
    let (a, b) = (
        settled(&fuller(), 1.0, 1.0),
        settled(&stress_response(), 1.0, 1.0),
    );
    assert!(level(&a, "p53") < level(&b, "p53"));
}

#[test]
fn p21_alone_does_not_stop_the_cycle_but_adds_to_hypoxia() {
    let stressed = settled(&fuller(), 1.0, 1.0);
    assert!(stressed.hypoxia_response() < QUIESCENCE_THRESHOLD);
    assert!(level(&stressed, "rb") > level(&settled(&fuller(), 1.0, 0.0), "rb"));
    // Just above the oxygen that arrests a calm cell, a stressed one stops.
    let o2 = arrest_oxygen(&fuller()) + 0.05;
    assert!(settled(&fuller(), o2, 0.0).hypoxia_response() < QUIESCENCE_THRESHOLD);
    assert!(settled(&fuller(), o2, 1.0).hypoxia_response() >= QUIESCENCE_THRESHOLD);
}

#[test]
fn the_cycle_settles_to_one_state_however_it_starts() {
    // Started arrested, RB fully on and E2F and cyclin E off, a cell ends
    // where one started cycling does: the switch has no memory, as HIF-1a
    // had none, so a cell that regains its oxygen divides again.
    let mut arrested = fuller();
    for node in &mut arrested.nodes {
        match node.name.as_str() {
            "rb" => node.initial = 1.0,
            "e2f" | "cyclin_e" => node.initial = 0.0,
            _ => {}
        }
    }
    for o2 in [1.0, 0.4, 0.32, 0.31, 0.0] {
        let (a, b) = (settled(&fuller(), o2, 0.0), settled(&arrested, o2, 0.0));
        assert!(
            (a.hypoxia_response() - b.hypoxia_response()).abs() < 1e-6,
            "at {o2}: RB {} against {}",
            a.hypoxia_response(),
            b.hypoxia_response()
        );
    }
}

#[test]
fn a_fuller_cell_starts_where_it_settles() {
    // The initial levels in the file are the well fed steady state, so the
    // first cell does not spend its first steps catching up.
    let spec = fuller();
    let net = settled(&spec, 1.0, 0.0);
    for (node, settled) in spec.nodes.iter().zip(&net.state().values) {
        if ["hif1", "vegf", "bcl2"].contains(&node.name.as_str()) {
            continue;
        }
        assert!(
            (node.initial - settled.0).abs() < 0.01,
            "{} starts at {} and settles at {}",
            node.name,
            node.initial,
            settled.0
        );
    }
}

#[test]
fn ras_stuck_on_raises_erk_as_far_through_raf_and_mek() {
    // The mutant RAS of chapter 2, laid over each network: switched on by
    // nothing rather than by what lies above it.
    let mutant: NetworkSpec = toml::from_str(include_str!(
        "../../../examples/networks/ras-erk-mutant.toml"
    ))
    .unwrap();
    let ras = mutant
        .nodes
        .iter()
        .find(|n| n.name == "ras")
        .unwrap()
        .clone();
    let stuck = |spec: NetworkSpec| {
        let mut spec = spec;
        for node in &mut spec.nodes {
            if node.name == "ras" {
                *node = ras.clone();
            }
        }
        level(&settled(&spec, 1.0, 0.0), "erk")
    };
    let (a, b) = (stuck(fuller()), stuck(simpler()));
    assert!((a - b).abs() < 0.01, "ERK {a:.4} against {b:.4}");
    assert!(a > 0.5, "well above the 0.33 of a normal RAS: {a}");
}
