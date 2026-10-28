//! A handful of numbers that say what shape a tissue has taken.
//!
//! A snapshot lists every cell, which is everything and says nothing at a
//! glance. What a hypothesis is judged on is a few numbers: is there a core,
//! how large, how deep does living tissue reach from its supply, how fast it
//! turns over. Measured here, once, so that two runs are compared on the same
//! geometry rather than each reader reimplementing it.

use cellweave_core::OxygenSource;

use crate::cpm::OBSTACLE;
use crate::simulation::{Fate, SimState, Simulation};

/// What a tissue looks like at one step.
///
/// The three radii are equivalent radii: the radius of a disc with the same
/// area. A tissue is never a perfect disc, and the equivalent radius is the
/// one number that does not depend on where its boundary happens to wobble.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Measurements {
    /// Index of the step last taken.
    pub mcs: u64,
    /// Cells alive and dividing.
    pub cycling: usize,
    /// Cells alive and arrested by hypoxia.
    pub quiescent: usize,
    /// Cells dead and not yet resorbed.
    pub necrotic: usize,
    /// Cells committed to differentiate and leaving the layer.
    pub differentiating: usize,
    /// Cells in apoptosis, breaking up.
    pub apoptotic: usize,
    /// Level of the lethal stress held at its source at this step, zero
    /// without an exposure or while it is off.
    pub exposure: f64,
    /// Cells ever born, the first one included, or the whole first sheet.
    pub born: u64,
    /// Cells that have left for good: resorbed to nothing after dying, or
    /// shed after differentiating.
    pub gone: u64,
    /// Clones with at least one living cell.
    pub clones: usize,
    /// Living cells descended from the scenario's lesion, zero without one.
    pub lesion: usize,
    /// Equivalent radius of everything the tissue covers, in pixels: living,
    /// necrotic and departing cells, and the vessel they surround if there is
    /// one.
    pub tissue_radius: f64,
    /// Equivalent radius of the necrotic area, in pixels. In a spheroid this
    /// is the core; around a vessel the dead lie outside the living, and it is
    /// only a measure of how much dead tissue there is.
    pub necrotic_radius: f64,
    /// How far living tissue reaches from where its oxygen comes from, in
    /// pixels.
    ///
    /// From a bath, oxygen arrives at the outer edge and the living tissue is
    /// a ring outside the core: the tissue radius less the core's. From a
    /// vessel, it arrives at the wall and the living tissue is a sleeve
    /// around it: the radius of vessel and living together, less the vessel's.
    /// Either way it is the depth diffusion can feed, which is the number a
    /// pathologist reads off a section. A sheet has no supply to measure
    /// from, and there the number means nothing.
    pub viable_depth: f64,
}

impl Measurements {
    /// Cells alive, cycling or quiescent.
    pub fn living(&self) -> usize {
        self.cycling + self.quiescent
    }
}

/// The radius of a disc of `area` pixels.
fn equivalent_radius(area: u64) -> f64 {
    (area as f64 / std::f64::consts::PI).sqrt()
}

impl<S: SimState> Simulation<S> {
    /// Measure the tissue as it stands.
    ///
    /// Reads the lattice and every cell's fate, and nothing else: it draws no
    /// random number and moves nothing, so measuring a run does not change it.
    pub fn measure(&self) -> Measurements {
        let lattice = self.lattice();
        let (mut cycling, mut quiescent, mut necrotic) = (0, 0, 0);
        let (mut differentiating, mut apoptotic) = (0, 0);
        let (mut living_area, mut necrotic_area, mut leaving_area) = (0u64, 0u64, 0u64);
        for (id, record) in lattice.cell_records() {
            if record.volume.0 == 0 {
                continue;
            }
            let area = u64::from(record.volume.0);
            match self.fate(id) {
                Some(Fate::Cycling) => {
                    cycling += 1;
                    living_area += area;
                }
                Some(Fate::Quiescent) => {
                    quiescent += 1;
                    living_area += area;
                }
                Some(Fate::Necrotic) => {
                    necrotic += 1;
                    necrotic_area += area;
                }
                Some(Fate::Differentiating) => {
                    differentiating += 1;
                    leaving_area += area;
                }
                Some(Fate::Apoptotic) => {
                    apoptotic += 1;
                    leaving_area += area;
                }
                None => {}
            }
        }

        let vessel_area = match self.scenario().oxygen_source {
            OxygenSource::Bath => 0,
            OxygenSource::Vessel { .. } => {
                let (w, h) = lattice.dims();
                (0..h)
                    .flat_map(|y| (0..w).map(move |x| (x, y)))
                    .filter(|&(x, y)| lattice.occupant(x, y) == OBSTACLE)
                    .count() as u64
            }
        };

        let tissue_radius =
            equivalent_radius(vessel_area + living_area + necrotic_area + leaving_area);
        let necrotic_radius = equivalent_radius(necrotic_area);
        let viable_depth = match self.scenario().oxygen_source {
            OxygenSource::Bath => tissue_radius - necrotic_radius,
            OxygenSource::Vessel { .. } => {
                equivalent_radius(vessel_area + living_area) - equivalent_radius(vessel_area)
            }
        };

        Measurements {
            mcs: self.last_mcs(),
            cycling,
            quiescent,
            necrotic,
            differentiating,
            apoptotic,
            exposure: self
                .scenario()
                .exposure
                .map_or(0.0, |e| e.level_at(self.last_mcs())),
            born: self.born(),
            gone: self.gone(),
            clones: self.clone_sizes().iter().filter(|&&n| n > 0).count(),
            lesion: self.lesion_cells(),
            tissue_radius,
            necrotic_radius,
            viable_depth,
        }
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]

    use cellweave_core::NetworkSpec;

    use super::*;
    use crate::simulation::{DEFAULT_SEED, Scenario};
    use cellweave_core::CellId;

    fn spec() -> NetworkSpec {
        toml::from_str(include_str!("../../../examples/networks/ras-erk.toml"))
            .expect("the shipped network parses")
    }

    fn scenario(oxygen_source: OxygenSource) -> Scenario {
        Scenario {
            width: 80,
            height: 80,
            temperature: 10.0,
            oxygen_source,
            oxygen_uptake: 0.0003,
            seed: DEFAULT_SEED,
            mutation_rate: 0.0,
            suppress: Vec::new(),
            rewiring: 0.0,
            seeding: Default::default(),
            exposure: None,
            lesion: None,
        }
    }

    #[test]
    fn a_lone_cell_in_a_bath_is_all_viable_depth() {
        let sim = Simulation::new(scenario(OxygenSource::Bath), &spec()).unwrap();
        let m = sim.measure();
        assert_eq!((m.cycling, m.quiescent, m.necrotic), (1, 0, 0));
        // Seeded as a disc of radius 8, so its equivalent radius is close to
        // that, and with no core all of it is viable.
        assert!((m.tissue_radius - 8.0).abs() < 0.5, "{}", m.tissue_radius);
        assert_eq!(m.necrotic_radius, 0.0);
        assert_eq!(m.viable_depth, m.tissue_radius);
    }

    #[test]
    fn around_a_vessel_the_depth_starts_at_the_wall() {
        let sim = Simulation::new(scenario(OxygenSource::Vessel { radius: 6 }), &spec()).unwrap();
        let m = sim.measure();
        let living = sim.lattice().volume(CellId(1)).unwrap().0;
        assert_eq!(m.living(), 1);
        // The tissue covers the vessel and the cell; the depth is what the
        // cell adds to the vessel's radius, so it is less than the radius the
        // cell would have on its own.
        assert!(m.tissue_radius > 6.0);
        assert!(m.viable_depth > 0.0);
        assert!(m.viable_depth < equivalent_radius(u64::from(living)));
    }

    #[test]
    fn measuring_does_not_change_the_run() {
        let run = |measure: bool| {
            let (mut sim, _) = Simulation::new(scenario(OxygenSource::Bath), &spec())
                .unwrap()
                .start()
                .unwrap();
            for _ in 0..40 {
                sim.step().unwrap();
                if measure {
                    sim.measure();
                }
            }
            format!("{:?}", sim.snapshot())
        };
        assert_eq!(run(true), run(false));
    }
}
