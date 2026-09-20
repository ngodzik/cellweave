//! Where the lattice and the field meet.
//!
//! Neither knows about the other. The lattice knows which cell holds each pixel
//! and what kind of cell it is; the field knows one concentration per pixel.
//! The three functions here carry information across, in both directions: what
//! the cells take up and where the bath can reach go down into the field, and
//! what each cell sits in comes back up, ready to be handed to its network.

use crate::cpm::CpmLattice;
use crate::pde::ScalarField;
use cellweave_core::{CellKind, DiffusionField};

/// Whether a cell of this kind draws on the field at all.
fn consumes(kind: CellKind) -> bool {
    !matches!(kind, CellKind::Medium | CellKind::Necrotic)
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
            let alive = lattice
                .cell_kind(cellweave_core::CellId(id))
                .is_some_and(consumes);
            if alive {
                field.set_uptake(x, y, rate);
            }
        }
    }
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
            let id = lattice.occupant(x, y) as usize;
            if id != 0 && id < n {
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
}
