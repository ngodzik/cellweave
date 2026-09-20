//! RAS-ERK-HIF-1α signaling network, integrated with RK4.
//!
//! Protein indices:
//!   0 RAS  (activated fraction, 0-1)
//!   1 ERK  (activated fraction, 0-1)
//!   2 HIF1 (hypoxia-inducible factor, 0-1)
//!   3 VEGF (secretion rate, 0-1)
//!   4 BCL2 (survival factor, 0-1)

use cellweave_core::{
    CellInputs, CellMechanics, Concentration, ProteinState, SignalingNetwork, TimeStep, Volume,
};

const N_PROTEINS: usize = 5;

/// Oxygen, as a fraction of the bath, below which a cell can no longer make
/// enough energy and starts taking damage.
///
/// Necrosis is not a decision the cell takes, it is a failure: below a critical
/// oxygen level oxidative metabolism stops, ATP falls, and membranes give way.
/// No survival programme rescues that, which is why the damage term below acts
/// on BCL2 directly and is not offset by HIF-1α.
const ANOXIA_THRESHOLD: f64 = 0.15;

/// How fast survival collapses at zero oxygen.
///
/// Set so that the steady state of BCL2 crosses the usual death threshold of
/// 0.25 around 5% of the bath, and stays intact above 10%: at 5%, anoxia is
/// 0.67 and 0.25 x 0.67 outweighs what ERK and HIF-1α feed in.
const ANOXIA_DAMAGE: f64 = 0.25;

// Protein indices.
const RAS: usize = 0;
const ERK: usize = 1;
const HIF1: usize = 2;
const VEGF_SEC: usize = 3;
const BCL2: usize = 4;

/// Minimal RAS-ERK-HIF-1α network.
///
/// Parameters are placeholder values; they will be calibrated from BioModels
/// SBML records in a later milestone.
#[derive(Debug, Clone)]
pub struct RasErkNetwork {
    /// Current protein state.
    state: ProteinState,
    /// Baseline target volume (modified by ERK).
    base_target_volume: u32,
}

impl RasErkNetwork {
    /// Create a network at steady state with no stimulation.
    pub fn new(base_target_volume: u32) -> Self {
        let mut state = ProteinState::zeros(N_PROTEINS);
        // Basal activities from literature estimates.
        state.values[RAS] = Concentration(0.05);
        state.values[ERK] = Concentration(0.05);
        state.values[HIF1] = Concentration(0.05);
        state.values[VEGF_SEC] = Concentration(0.05);
        state.values[BCL2] = Concentration(0.70);
        Self {
            state,
            base_target_volume,
        }
    }

    /// Reference volume of this lineage, which the target volume and the
    /// division threshold are both measured against.
    ///
    /// A daughter inherits it unchanged, so the size at which cells divide is a
    /// property of the lineage rather than of each cell's own birth size, which
    /// is how a cell type keeps a characteristic size across generations.
    pub fn base_target_volume(&self) -> u32 {
        self.base_target_volume
    }

    /// The survival signal, BCL2. A cell whose survival collapses is necrotic.
    pub fn survival(&self) -> f64 {
        self.state.values[BCL2].0
    }

    /// The hypoxia response, HIF-1α. It rises as oxygen falls, and a cell whose
    /// response is high has stopped cycling: hypoxia arrests the cell cycle
    /// well before it kills.
    pub fn hypoxia_response(&self) -> f64 {
        self.state.values[HIF1].0
    }

    /// ODE right-hand side: dX/dt = f(X, inputs).
    fn derivatives(x: &[f64; N_PROTEINS], inputs: &CellInputs) -> [f64; N_PROTEINS] {
        let [ras, erk, hif1, vegf, bcl2] = *x;

        // EGF → RAS (Michaelis-Menten activation, 1st-order degradation).
        let egf = inputs.egf.0;
        let dras = 0.5 * egf * (1.0 - ras) - 0.3 * ras;

        // RAS → ERK cascade.
        let derk = 0.8 * ras * (1.0 - erk) - 0.4 * erk;

        // Hypoxia → HIF-1α (O2 suppresses HIF-1α).
        let o2 = inputs.o2.0.clamp(0.0, 1.0);
        let dhif1 = 0.6 * (1.0 - o2) * (1.0 - hif1) - 0.5 * hif1;

        // HIF-1α → VEGF secretion.
        let dvegf = 0.7 * hif1 - 0.3 * vegf;

        // ERK ↑ BCL2, HIF-1α ↑ BCL2 (survival under stress), and below the
        // anoxia threshold damage pulls it down regardless of either.
        let anoxia = ((ANOXIA_THRESHOLD - o2) / ANOXIA_THRESHOLD).max(0.0);
        let dbcl2 = 0.2 * erk + 0.1 * hif1 - 0.15 * (bcl2 - 0.5) - ANOXIA_DAMAGE * anoxia;

        [dras, derk, dhif1, dvegf, dbcl2]
    }

    /// Classic RK4 step over `dt` seconds.
    fn rk4(&mut self, dt: f64, inputs: &CellInputs) {
        let x: [f64; N_PROTEINS] = std::array::from_fn(|i| self.state.values[i].0);

        let k1 = Self::derivatives(&x, inputs);
        let x2: [f64; N_PROTEINS] = std::array::from_fn(|i| x[i] + 0.5 * dt * k1[i]);
        let k2 = Self::derivatives(&x2, inputs);
        let x3: [f64; N_PROTEINS] = std::array::from_fn(|i| x[i] + 0.5 * dt * k2[i]);
        let k3 = Self::derivatives(&x3, inputs);
        let x4: [f64; N_PROTEINS] = std::array::from_fn(|i| x[i] + dt * k3[i]);
        let k4 = Self::derivatives(&x4, inputs);

        for i in 0..N_PROTEINS {
            let new_val = x[i] + (dt / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
            self.state.values[i] = Concentration(new_val.clamp(0.0, 1.0));
        }
    }
}

impl SignalingNetwork for RasErkNetwork {
    fn step(&mut self, dt: TimeStep, inputs: &CellInputs) {
        self.rk4(dt.0, inputs);
    }

    fn state(&self) -> &ProteinState {
        &self.state
    }

    fn mechanics(&self) -> CellMechanics {
        let erk = self.state.values[ERK].0;
        // High ERK pushes the cell toward proliferation. The coefficient is set so
        // that a well stimulated cell aims at roughly twice its newborn size, which
        // is what lets it reach a division threshold at all: ERK saturates near
        // 0.55 in this network, so a coefficient of 0.5 would cap growth below any
        // sensible threshold and no cell would ever divide.
        let target = (self.base_target_volume as f64 * (1.0 + 1.5 * erk)) as u32;
        // Low BCL2 → weaker survival → let volume shrink (apoptosis onset).
        let _bcl2 = self.state.values[BCL2].0;
        CellMechanics {
            target_volume: Volume(target),
            lambda_volume: 50.0,
            j_medium: 16.0,
            // Low E-cad surrogate: cells with high ERK become more invasive.
            j_self: 2.0 + 8.0 * (1.0 - erk),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn zero_inputs() -> CellInputs {
        CellInputs {
            o2: Concentration(1.0),
            egf: Concentration(0.0),
            vegf: Concentration(0.0),
            drugs: vec![],
        }
    }

    #[test]
    fn no_stimulation_stays_bounded() {
        let mut net = RasErkNetwork::new(100);
        let inputs = zero_inputs();
        for _ in 0..1000 {
            net.step(TimeStep(0.01), &inputs);
        }
        for c in &net.state().values {
            assert!((0.0..=1.0).contains(&c.0), "out of bounds: {}", c.0);
        }
    }

    #[test]
    fn egf_elevates_erk() {
        let mut net = RasErkNetwork::new(100);
        let inputs = CellInputs {
            o2: Concentration(1.0),
            egf: Concentration(1.0),
            vegf: Concentration(0.0),
            drugs: vec![],
        };
        for _ in 0..500 {
            net.step(TimeStep(0.05), &inputs);
        }
        let erk = net.state().values[ERK].0;
        assert!(erk > 0.5, "ERK should be elevated with EGF, got {:.3}", erk);
    }

    #[test]
    fn hypoxia_elevates_hif1() {
        let mut net = RasErkNetwork::new(100);
        let inputs = CellInputs {
            o2: Concentration(0.0),
            egf: Concentration(0.0),
            vegf: Concentration(0.0),
            drugs: vec![],
        };
        for _ in 0..500 {
            net.step(TimeStep(0.05), &inputs);
        }
        let hif1 = net.state().values[HIF1].0;
        // Steady-state for dHIF = 0.6*(1-hif) - 0.5*hif → hif = 0.6/1.1 ≈ 0.545.
        assert!(
            hif1 > 0.5,
            "HIF-1α should be elevated under hypoxia, got {:.3}",
            hif1
        );
    }

    /// The binary's inputs, at a given oxygen level.
    fn at_oxygen(o2: f64) -> CellInputs {
        CellInputs {
            o2: Concentration(o2),
            egf: Concentration(0.2),
            vegf: Concentration(0.0),
            drugs: vec![],
        }
    }

    fn settled_at(o2: f64) -> RasErkNetwork {
        let mut net = RasErkNetwork::new(200);
        for _ in 0..400 {
            net.step(TimeStep(1.0), &at_oxygen(o2));
        }
        net
    }

    #[test]
    fn sustained_anoxia_collapses_survival() {
        // Below the anoxia threshold damage outweighs everything ERK and HIF-1α
        // feed into survival. At 3% of the bath the cell is dead outright.
        assert!(settled_at(0.03).survival() < 0.05);
        // At 5% it has crossed the usual death threshold of 0.25.
        assert!(settled_at(0.05).survival() < 0.25);
    }

    #[test]
    fn hypoxia_short_of_anoxia_leaves_survival_intact() {
        // Hypoxic, with a strong HIF-1α response, but still making energy: the
        // cell stops cycling, it does not die.
        let net = settled_at(0.15);
        assert!(net.survival() > 0.9, "survival {:.3}", net.survival());
        assert!(net.hypoxia_response() > 0.4);
    }

    #[test]
    fn the_hypoxia_response_rises_as_oxygen_falls() {
        let levels = [1.0, 0.5, 0.2, 0.05];
        let responses: Vec<f64> = levels
            .iter()
            .map(|&o2| settled_at(o2).hypoxia_response())
            .collect();
        assert!(
            responses.windows(2).all(|w| w[0] < w[1]),
            "not monotone: {responses:?}"
        );
        assert!(responses[0] < 0.05, "no response under full oxygen");
    }

    #[test]
    fn the_network_answers_well_within_a_cell_cycle() {
        // Three clocks have to agree: the field is relaxed to equilibrium every
        // lattice step, and a cell grows to division in about twenty lattice
        // steps. The network settles in about fifteen units of its own time, so
        // the binary gives it eight units per lattice step and it answers in two.
        // This pins the fifteen: were it to drift toward the cell cycle, a cell's
        // fate would trail its surroundings by a generation.
        for o2 in [0.8, 0.02] {
            let target = settled_at(o2);
            let mut net = RasErkNetwork::new(200);
            let mut settled_after = None;
            for step in 1..=100 {
                net.step(TimeStep(1.0), &at_oxygen(o2));
                let worst = (0..N_PROTEINS)
                    .map(|i| (net.state().values[i].0 - target.state().values[i].0).abs())
                    .fold(0.0_f64, f64::max);
                if worst < 0.01 {
                    settled_after = Some(step);
                    break;
                }
            }
            let steps = settled_after.expect("the network never settled");
            assert!(steps < 40, "settled after {steps} steps at oxygen {o2}");
        }
    }
}
