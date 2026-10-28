//! Where the lattice and the field meet.
//!
//! Neither knows about the other. The lattice knows which cell holds each pixel
//! and what kind of cell it is; the field knows one concentration per pixel.
//! The three functions here carry information across, in both directions: what
//! the cells take up and where the bath can reach go down into the field, and
//! what each cell sits in comes back up, ready to be handed to its network.

use crate::cpm::{CpmLattice, OBSTACLE};
use crate::pde::ScalarField;
use cellweave_core::{CellKind, DiffusionField};

/// Whether a cell of this kind draws on the field at all. A dead cell does
/// not, and neither does one breaking up in apoptosis.
fn consumes(kind: CellKind) -> bool {
    !matches!(
        kind,
        CellKind::Medium | CellKind::Necrotic | CellKind::Apoptotic
    )
}

/// Make every pixel of every living cell a sink of `rate`.
///
/// The field's sinks are rewritten from scratch, so a cell that died or moved
/// since the last call stops consuming where it no longer is.
pub fn apply_uptake(field: &mut ScalarField, lattice: &CpmLattice, rate: f64) {
    field.clear_uptake();
    let (w, h) = lattice.dims();
    for y in 0..h {
        for x in 0..w {
            let id = lattice.occupant(x, y);
            if id == 0 || id == OBSTACLE {
                continue;
            }
            let alive = lattice
                .cell_kind(cellweave_core::CellId(id))
                .is_some_and(consumes);
            if alive {
                field.set_uptake(x, y, rate);
            }
        }
    }
}

/// Make every medium pixel a sink of `rate` as well.
///
/// Around a vessel there is no empty liquid; there is stroma, ordinary tissue
/// that consumes too. Without this the medium beyond the tumour would fill up
/// from the vessel and hand oxygen back to cells that should be starving. Call
/// after [`apply_uptake`], which clears the sinks first.
pub fn apply_background_uptake(field: &mut ScalarField, lattice: &CpmLattice, rate: f64) {
    let (w, h) = lattice.dims();
    for y in 0..h {
        for x in 0..w {
            if lattice.occupant(x, y) == 0 {
                field.set_uptake(x, y, rate);
            }
        }
    }
}

/// Make every pixel of every cell secrete what `per_cell` says that cell
/// does, indexed by identifier.
///
/// The field's sources are rewritten from scratch, like its sinks, so a cell
/// that died or moved since the last call stops secreting where it no longer
/// is. A cell missing from `per_cell` secretes nothing.
pub fn apply_production(field: &mut ScalarField, lattice: &CpmLattice, per_cell: &[f64]) {
    field.clear_production();
    let (w, h) = lattice.dims();
    for y in 0..h {
        for x in 0..w {
            let id = lattice.occupant(x, y);
            if id == 0 || id == OBSTACLE {
                continue;
            }
            if let Some(&rate) = per_cell.get(id as usize)
                && rate > 0.0
            {
                field.set_production(x, y, rate);
            }
        }
    }
}

/// Hold every obstacle pixel at `value`, which makes the obstacle a source.
///
/// A vessel is an obstacle that emits: cells cannot enter it, and it keeps the
/// oxygen at its wall at a fixed level. This is the geometry of a tumour cord,
/// tissue living around a capillary out to the depth diffusion can feed.
pub fn hold_obstacles_at(field: &mut ScalarField, lattice: &CpmLattice, value: f64) {
    field.hold_pixels(value, |x, y| lattice.occupant(x, y) == OBSTACLE);
}

/// Hold every obstacle pixel `carries` says carries blood at `value`.
///
/// A vessel that has grown is held only where blood flows through it; a
/// blind sprout is an obstacle that supplies nothing. With every vessel
/// carrying blood this is [`hold_obstacles_at`].
pub fn hold_vessels_at(
    field: &mut ScalarField,
    lattice: &CpmLattice,
    value: f64,
    carries: impl Fn(u32, u32) -> bool,
) {
    field.hold_pixels(value, |x, y| {
        lattice.occupant(x, y) == OBSTACLE && carries(x, y)
    });
}

/// Hold every obstacle pixel at `value` times what `delivery` says it
/// delivers there, between zero and one, where that is more than the field
/// already holds around it; leave it free otherwise.
///
/// A vessel delivers in proportion to the blood running through it, so a
/// vessel network is a source whose strength varies along it. It gives what
/// its blood brings and takes nothing: held below the tissue around it, a
/// vessel carrying a trickle would draw the oxygen out of the tissue as if
/// its blood could carry any amount away. So where the pixels around a
/// vessel pixel, as the field last stood, already hold more than it would,
/// it is left free, as a vessel carrying no blood is. With every vessel
/// delivering in full, and nothing richer than a vessel, this is
/// [`hold_obstacles_at`].
pub fn hold_vessels_by(
    field: &mut ScalarField,
    lattice: &CpmLattice,
    value: f64,
    delivery: impl Fn(u32, u32) -> f64,
) {
    let (w, h) = lattice.dims();
    let mut holds = vec![None; w as usize * h as usize];
    for y in 0..h {
        for x in 0..w {
            if lattice.occupant(x, y) != OBSTACLE {
                continue;
            }
            let share = delivery(x, y);
            if share <= 0.0 {
                continue;
            }
            let held = value * share.min(1.0);
            let around: Vec<f64> = [
                (x.wrapping_sub(1), y),
                (x + 1, y),
                (x, y.wrapping_sub(1)),
                (x, y + 1),
            ]
            .iter()
            .filter(|&&(i, j)| i < w && j < h && lattice.occupant(i, j) != OBSTACLE)
            .map(|&(i, j)| field.concentration_at(i, j))
            .collect();
            let richer =
                !around.is_empty() && around.iter().sum::<f64>() / around.len() as f64 >= held;
            if !richer {
                holds[y as usize * w as usize + x as usize] = Some(held);
            }
        }
    }
    field.hold_pixels_by(|x, y| holds[y as usize * w as usize + x as usize]);
}

/// Hold the medium the bath can reach at `value`.
///
/// Reachability is decided by [`ScalarField::hold_bath`]; this only says which
/// pixels are open, namely those no cell occupies.
pub fn hold_medium_bath(field: &mut ScalarField, lattice: &CpmLattice, value: f64) {
    field.hold_bath(value, |x, y| lattice.occupant(x, y) == 0);
}

/// Mean concentration over each cell's pixels, indexed by identifier.
///
/// Index 0 is the medium and is left at zero. A cell with no pixels, which
/// bookkeeping should never produce, reads as zero rather than dividing by it.
pub fn mean_per_cell(field: &ScalarField, lattice: &CpmLattice) -> Vec<f64> {
    let n = lattice.cell_count() + 1;
    let mut sum = vec![0.0; n];
    let mut count = vec![0u32; n];
    let (w, h) = lattice.dims();
    for y in 0..h {
        for x in 0..w {
            let id = lattice.occupant(x, y);
            if id == 0 || id == OBSTACLE {
                continue;
            }
            let id = id as usize;
            if id < n {
                sum[id] += field.concentration_at(x, y);
                count[id] += 1;
            }
        }
    }
    sum.iter()
        .zip(&count)
        .map(|(s, c)| if *c == 0 { 0.0 } else { s / f64::from(*c) })
        .collect()
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]

    use super::*;
    use cellweave_core::{CellId, Pos2, TimeStep};

    /// A packed cluster: one cell in the middle, ringed by others.
    fn cluster() -> (CpmLattice, CellId, Vec<CellId>) {
        let mut lat = CpmLattice::new(80, 80, 10.0).unwrap();
        let centre = lat.add_cell(CellKind::Tumor, Pos2 { x: 40, y: 40 }, 6);
        let ring: Vec<CellId> = (0..8)
            .map(|k| {
                let a = f64::from(k) * std::f64::consts::TAU / 8.0;
                let x = (40.0 + 13.0 * a.cos()).round() as u32;
                let y = (40.0 + 13.0 * a.sin()).round() as u32;
                lat.add_cell(CellKind::Tumor, Pos2 { x, y }, 6)
            })
            .collect();
        (lat, centre, ring)
    }

    fn settled_field(lat: &CpmLattice) -> ScalarField {
        let mut field = ScalarField::new(80, 80, 0.2, 0.0, 1.0).unwrap();
        apply_uptake(&mut field, lat, 0.003);
        hold_medium_bath(&mut field, lat, 1.0);
        field.relax(TimeStep(1.0), 1e-7, 50_000).unwrap();
        field
    }

    #[test]
    fn a_cell_in_the_middle_of_a_cluster_reads_less_oxygen_than_one_on_the_rim() {
        // No network involved: this is the lattice and the field alone, which is
        // what makes the reading a property of geometry rather than of biology.
        let (lat, centre, ring) = cluster();
        let field = settled_field(&lat);

        let o2 = mean_per_cell(&field, &lat);
        let inner = o2[centre.0 as usize];
        let outer = ring.iter().map(|id| o2[id.0 as usize]).fold(0.0, f64::max);

        assert!(
            inner < outer,
            "centre {inner:.3} should be below rim {outer:.3}"
        );
        assert!(
            inner < 0.9 && outer > inner + 0.05,
            "gradient too weak to mean anything"
        );
    }

    #[test]
    fn a_necrotic_cell_stops_consuming() {
        let (mut lat, centre, _) = cluster();
        let starved = mean_per_cell(&settled_field(&lat), &lat)[centre.0 as usize];

        lat.set_kind(centre, CellKind::Necrotic);
        let relieved = mean_per_cell(&settled_field(&lat), &lat)[centre.0 as usize];

        assert!(
            relieved > starved,
            "a dead centre no longer draws: {starved:.3} then {relieved:.3}"
        );
    }

    #[test]
    fn the_medium_reads_zero_and_every_cell_reads_a_mean() {
        let (lat, _, _) = cluster();
        let o2 = mean_per_cell(&settled_field(&lat), &lat);

        assert_eq!(o2.len(), lat.cell_count() + 1);
        assert_eq!(o2[0], 0.0);
        assert!(o2[1..].iter().all(|c| *c > 0.0 && *c <= 1.0));
    }

    #[test]
    fn a_vessel_feeds_the_cell_beside_it_and_not_the_one_far_away() {
        let mut lat = CpmLattice::new(80, 80, 10.0).unwrap();
        lat.add_obstacle_disc(Pos2 { x: 40, y: 40 }, 4);
        let near = lat.add_cell(CellKind::Tumor, Pos2 { x: 50, y: 40 }, 5);
        let far = lat.add_cell(CellKind::Tumor, Pos2 { x: 70, y: 40 }, 5);

        let mut field = ScalarField::new(80, 80, 0.2, 0.0, 0.0).unwrap();
        apply_uptake(&mut field, &lat, 0.003);
        hold_obstacles_at(&mut field, &lat, 1.0);
        field.relax(TimeStep(1.0), 1e-7, 50_000).unwrap();

        let o2 = mean_per_cell(&field, &lat);
        let (n, f) = (o2[near.0 as usize], o2[far.0 as usize]);
        assert!(
            n > f + 0.1,
            "near {n:.3} should read far more than far {f:.3}"
        );
        // The obstacle itself has no entry to read.
        assert_eq!(o2.len(), lat.cell_count() + 1);
    }

    #[test]
    fn a_vessel_delivers_what_its_share_says_and_a_dry_one_nothing() {
        let mut lat = CpmLattice::new(10, 1, 10.0).unwrap();
        for x in 0..3 {
            lat.add_obstacle_disc(Pos2 { x, y: 0 }, 0);
        }
        let mut field = ScalarField::new(10, 1, 0.2, 0.0, 0.0).unwrap();
        hold_vessels_by(&mut field, &lat, 2.0, |x, _| [0.0, 0.5, 3.0][x as usize]);
        assert_eq!(
            field.concentration_at(0, 0),
            0.0,
            "a dry vessel is not held"
        );
        assert_eq!(field.concentration_at(1, 0), 1.0, "half of 2");
        assert_eq!(field.concentration_at(2, 0), 2.0, "never more than in full");
    }

    #[test]
    fn a_vessel_poorer_than_the_tissue_around_it_takes_nothing_from_it() {
        let mut lat = CpmLattice::new(10, 1, 10.0).unwrap();
        lat.add_obstacle_disc(Pos2 { x: 5, y: 0 }, 0);
        let mut field = ScalarField::new(10, 1, 0.2, 0.0, 0.0).unwrap();
        field.fill(0.9);
        hold_vessels_by(&mut field, &lat, 1.0, |_, _| 0.5);
        assert_eq!(
            field.concentration_at(5, 0),
            0.9,
            "left free, not held down to half"
        );
        hold_vessels_by(&mut field, &lat, 1.0, |_, _| 1.0);
        assert_eq!(field.concentration_at(5, 0), 1.0, "held where it is richer");
    }
}
