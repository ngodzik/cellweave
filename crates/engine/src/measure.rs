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
    /// Pixels of vessel: the one the scenario laid, and whatever has grown
    /// from it.
    pub vessel: u64,
    /// Pixels of vessel that carry blood.
    pub perfused: u64,
    /// Sprouts still growing.
    pub sprouts: usize,
    /// Sprouts that have regressed so far, for carrying no blood too long.
    pub regressed: usize,
    /// Immune cells in the tissue, not counting those dying.
    pub immune: usize,
    /// Tumour cells immune cells have killed so far.
    pub killed: u64,
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
    /// From a capillary lying across the box, it is a layer along one side of
    /// it: the living area over the number of columns living cells hold.
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

/// Share of the living pixels [`Simulation::depth_from_vessels`] is read at
/// when a run reports it: nearly all of them, but not the odd cell squeezed
/// out to where no vessel feeds it.
pub const DEPTH_SHARE: f64 = 0.95;

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
        let mut living = vec![false; lattice.cell_count() + 1];
        for (id, record) in lattice.cell_records() {
            if record.volume.0 == 0 {
                continue;
            }
            let area = u64::from(record.volume.0);
            match self.fate(id) {
                Some(Fate::Cycling) => {
                    cycling += 1;
                    living_area += area;
                    living[id.0 as usize] = true;
                }
                Some(Fate::Quiescent) => {
                    quiescent += 1;
                    living_area += area;
                    living[id.0 as usize] = true;
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

        // A bath has no vessel, and nothing grows one there.
        let (vessel_area, perfused) = match self.scenario().oxygen_source {
            OxygenSource::Bath => (0, 0),
            OxygenSource::Vessel { .. } | OxygenSource::Capillary { .. } => {
                let (w, h) = lattice.dims();
                let vessel: Vec<(u32, u32)> = (0..h)
                    .flat_map(|y| (0..w).map(move |x| (x, y)))
                    .filter(|&(x, y)| lattice.occupant(x, y) == OBSTACLE)
                    .collect();
                let perfused = vessel
                    .iter()
                    .filter(|&&(x, y)| self.carries_blood(x, y))
                    .count();
                (vessel.len() as u64, perfused as u64)
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
            OxygenSource::Capillary { .. } => {
                let (w, h) = lattice.dims();
                let columns = (0..w)
                    .filter(|&x| {
                        (0..h).any(|y| {
                            let id = lattice.occupant(x, y);
                            id != OBSTACLE && living.get(id as usize).copied().unwrap_or(false)
                        })
                    })
                    .count();
                if columns == 0 {
                    0.0
                } else {
                    living_area as f64 / columns as f64
                }
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
            vessel: vessel_area,
            perfused,
            sprouts: self.sprouts().len(),
            regressed: self.regressed(),
            immune: self.immune_cells(),
            killed: self.killed(),
            tissue_radius,
            necrotic_radius,
            viable_depth,
        }
    }
}

impl<S: SimState> Simulation<S> {
    /// How far living tissue reaches from its nearest vessel, in pixels: the
    /// distance within which `share` of the living pixels lie. `None` when
    /// there is no vessel or nothing living.
    ///
    /// A tumour cord is a sleeve of living tissue around a vessel, as thick
    /// as diffusion can feed. Where vessels grow, this says whether the
    /// tissue grew thicker sleeves or more of them.
    pub fn depth_from_vessels(&self, share: f64) -> Option<f64> {
        let lattice = self.lattice();
        let (w, h) = lattice.dims();
        let pixels = || (0..h).flat_map(move |y| (0..w).map(move |x| (x, y)));
        let vessel: Vec<(i64, i64)> = pixels()
            .filter(|&(x, y)| lattice.occupant(x, y) == OBSTACLE)
            .map(|(x, y)| (i64::from(x), i64::from(y)))
            .collect();
        if vessel.is_empty() {
            return None;
        }
        let mut depths: Vec<f64> = pixels()
            .filter(|&(x, y)| {
                let id = lattice.occupant(x, y);
                id != 0
                    && id != OBSTACLE
                    && matches!(
                        self.fate(cellweave_core::CellId(id)),
                        Some(Fate::Cycling | Fate::Quiescent)
                    )
            })
            .map(|(x, y)| {
                let (x, y) = (i64::from(x), i64::from(y));
                let nearest = vessel
                    .iter()
                    .map(|&(vx, vy)| (vx - x).pow(2) + (vy - y).pow(2))
                    .min()
                    .unwrap_or(0);
                (nearest as f64).sqrt()
            })
            .collect();
        if depths.is_empty() {
            return None;
        }
        depths.sort_by(f64::total_cmp);
        let at = ((depths.len() - 1) as f64 * share.clamp(0.0, 1.0)).round() as usize;
        Some(depths[at])
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
            suppress_before: None,
            rewiring: 0.0,
            seeding: Default::default(),
            exposure: None,
            lesion: None,
            revert: None,
            angiogenesis: None,
            immunity: None,
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
