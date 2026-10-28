//! Vessels that grow towards the tissue that calls for them.
//!
//! A vessel here is a set of obstacle pixels held at full oxygen. Cells
//! starved of oxygen secrete VEGF, which diffuses and decays; where it is
//! high enough at a vessel's wall, the wall sends out a sprout. The sprout's
//! tip climbs the VEGF gradient a fraction of a pixel per step and lays a
//! disc of new vessel each time it reaches a new pixel, carving it out of
//! whatever cells were there. A tip stops when the call it follows falls
//! below the threshold, when it reaches the edge of the box, or when it runs
//! into another vessel and joins it.
//!
//! Each sprout remembers what it grew from. Blood flows through a new vessel
//! either at once or, when the scenario says so, only once its tip has
//! joined another vessel: then the sprout, and every sprout it grew from,
//! carries blood, and so does the one it joined, and everything that one
//! grew from.
//!
//! Nothing here is drawn at random. Which wall sprouts is the one reading the
//! most VEGF, and where a tip goes is where the gradient points, so a run
//! with vessels that may grow is the run without them, to the last bit, until
//! the first sprout lays its first pixel.

use cellweave_core::{Angiogenesis, DiffusionField, Perfusion, Pos2};

use crate::cpm::{CpmLattice, OBSTACLE};
use crate::pde::ScalarField;

/// Radius of the vessel a tip lays behind it, in pixels.
///
/// A capillary is five to ten micrometres across and a pixel is about one,
/// so a disc of radius two, five pixels across, is the narrowest vessel there
/// is, which is what a new sprout is.
pub const VESSEL_RADIUS: u32 = 2;

/// Closest a new sprout may start to a tip still growing, or to where
/// another sprout once started, in pixels: about a cell's diameter.
///
/// It stands in for the way tip cells keep the endothelial cells next to
/// them from becoming tips themselves. Without it every wall pixel reading
/// enough VEGF would sprout at once, into a solid front of vessel.
const SPROUT_SPACING: f64 = 16.0;

/// How many of the discs a tip laid last count as its own trail when it
/// looks ahead for another vessel. A tip turns slowly, so older discs are
/// out of its reach.
const TRAIL_KEPT: usize = 4;

/// A tip only looks for another vessel once this far from where it
/// sprouted, in pixels, clear of the vessel it came from.
const CLEAR_OF_ORIGIN: f64 = 2.0 * VESSEL_RADIUS as f64 + 2.0;

/// The leading end of a sprout.
#[derive(Debug, Clone, PartialEq)]
pub struct Tip {
    /// Where it is, in pixels, along the width of the box.
    pub x: f64,
    /// Where it is, in pixels, along the height of the box.
    pub y: f64,
    /// The wall pixel it sprouted from.
    origin: (u32, u32),
    /// Centres of the last few discs of vessel it laid, the latest last, or
    /// its origin before it has laid any.
    trail: Vec<(u32, u32)>,
    /// The sprout it leads.
    sprout: usize,
}

impl Tip {
    fn new(wall: (u32, u32), sprout: usize) -> Self {
        Self {
            x: f64::from(wall.0),
            y: f64::from(wall.1),
            origin: wall,
            trail: vec![wall],
            sprout,
        }
    }
}

/// The vessel one tip laid.
#[derive(Debug, Clone)]
struct Sprout {
    /// The sprout it grew from, or `None` for the vessel the scenario laid.
    parent: Option<usize>,
    /// Whether blood flows through it.
    perfused: bool,
}

/// Which sprout laid each pixel of vessel, and which sprouts carry blood.
#[derive(Debug, Clone, Default)]
struct Tree {
    sprouts: Vec<Sprout>,
    /// By pixel, row by row: the sprout that laid it, or `None` for the
    /// vessel the scenario laid and for every pixel that is not vessel.
    /// Empty until the first sprout.
    owner: Vec<Option<usize>>,
    width: u32,
}

impl Tree {
    fn index(&self, (x, y): (u32, u32)) -> usize {
        y as usize * self.width as usize + x as usize
    }

    fn owner(&self, pixel: (u32, u32)) -> Option<usize> {
        self.owner.get(self.index(pixel)).copied().flatten()
    }

    /// Lay a disc of vessel for `sprout` at `centre`, and remember which
    /// pixels it made vessel.
    fn lay(&mut self, lattice: &mut CpmLattice, centre: (u32, u32), sprout: usize) {
        let (w, h) = lattice.dims();
        let r = VESSEL_RADIUS;
        let r2 = i64::from(r * r);
        let rows = centre.1.saturating_sub(r)..(centre.1 + r + 1).min(h);
        let cols = centre.0.saturating_sub(r)..(centre.0 + r + 1).min(w);
        let fresh: Vec<(u32, u32)> = rows
            .flat_map(|y| cols.clone().map(move |x| (x, y)))
            .filter(|&(x, y)| {
                let (dx, dy) = (
                    i64::from(x) - i64::from(centre.0),
                    i64::from(y) - i64::from(centre.1),
                );
                dx * dx + dy * dy <= r2 && lattice.occupant(x, y) != OBSTACLE
            })
            .collect();
        lattice.add_obstacle_disc(
            Pos2 {
                x: centre.0,
                y: centre.1,
            },
            r,
        );
        for pixel in fresh {
            let i = self.index(pixel);
            self.owner[i] = Some(sprout);
        }
    }

    /// Blood through `sprout` and every sprout it grew from.
    fn perfuse(&mut self, sprout: usize) {
        let mut at = Some(sprout);
        while let Some(i) = at {
            self.sprouts[i].perfused = true;
            at = self.sprouts[i].parent;
        }
    }

    fn carries_blood(&self, pixel: (u32, u32)) -> bool {
        self.owner(pixel).is_none_or(|s| self.sprouts[s].perfused)
    }
}

/// Every sprout still growing, where every sprout ever started, and which of
/// them carry blood.
#[derive(Debug, Clone, Default)]
pub(crate) struct Vessels {
    tips: Vec<Tip>,
    origins: Vec<(u32, u32)>,
    tree: Tree,
}

impl Vessels {
    /// Whether the vessel at `(x, y)` carries blood: the vessel the scenario
    /// laid always does, a sprout once it is perfused. Meaningless for a
    /// pixel that is not vessel.
    pub(crate) fn carries_blood(&self, x: u32, y: u32) -> bool {
        self.tree.carries_blood((x, y))
    }

    /// The tips still growing.
    pub(crate) fn tips(&self) -> &[Tip] {
        &self.tips
    }

    /// How many sprouts have started, growing or not.
    pub(crate) fn sprouted(&self) -> usize {
        self.origins.len()
    }

    /// Move every tip, then let the wall reading the most VEGF sprout, at
    /// most once. Returns whether any vessel was laid, which is when cells
    /// may have lost pixels to it.
    pub(crate) fn grow(
        &mut self,
        lattice: &mut CpmLattice,
        vegf: &ScalarField,
        angiogenesis: &Angiogenesis,
    ) -> bool {
        if self.tree.owner.is_empty() {
            let (w, h) = lattice.dims();
            self.tree.owner = vec![None; w as usize * h as usize];
            self.tree.width = w;
        }
        let mut laid = false;
        let tree = &mut self.tree;
        self.tips
            .retain_mut(|tip| advance(tip, lattice, vegf, angiogenesis, tree, &mut laid));
        if let Some(wall) = self.strongest_free_wall(lattice, vegf, angiogenesis.threshold) {
            let sprout = self.tree.sprouts.len();
            self.tree.sprouts.push(Sprout {
                parent: self.tree.owner(wall),
                perfused: angiogenesis.perfusion == Perfusion::AtOnce,
            });
            self.tips.push(Tip::new(wall, sprout));
            self.origins.push(wall);
        }
        laid
    }

    /// The vessel pixel open to the tissue that reads the most VEGF, at
    /// least `threshold`, and is far enough from every tip and every origin.
    /// Ties go to the first in reading order.
    fn strongest_free_wall(
        &self,
        lattice: &CpmLattice,
        vegf: &ScalarField,
        threshold: f64,
    ) -> Option<(u32, u32)> {
        let (w, h) = lattice.dims();
        let near = |x: u32, y: u32, (ax, ay): (f64, f64)| {
            (f64::from(x) - ax).hypot(f64::from(y) - ay) < SPROUT_SPACING
        };
        let mut best: Option<((u32, u32), f64)> = None;
        for y in 1..h.saturating_sub(1) {
            for x in 1..w.saturating_sub(1) {
                if lattice.occupant(x, y) != OBSTACLE {
                    continue;
                }
                let level = vegf.concentration_at(x, y);
                if level < threshold || best.is_some_and(|(_, b)| level <= b) {
                    continue;
                }
                let open = [(x - 1, y), (x + 1, y), (x, y - 1), (x, y + 1)]
                    .iter()
                    .any(|&(i, j)| lattice.occupant(i, j) != OBSTACLE);
                let crowded = self.tips.iter().any(|t| near(x, y, (t.x, t.y)))
                    || self
                        .origins
                        .iter()
                        .any(|&(ox, oy)| near(x, y, (f64::from(ox), f64::from(oy))));
                if open && !crowded {
                    best = Some(((x, y), level));
                }
            }
        }
        best.map(|(wall, _)| wall)
    }
}

/// Move one tip up the gradient, laying vessel when it reaches a new pixel.
/// Returns whether it is still growing.
fn advance(
    tip: &mut Tip,
    lattice: &mut CpmLattice,
    vegf: &ScalarField,
    angiogenesis: &Angiogenesis,
    tree: &mut Tree,
    laid: &mut bool,
) -> bool {
    let (w, h) = lattice.dims();
    let (ix, iy) = (tip.x.round() as u32, tip.y.round() as u32);
    if vegf.concentration_at(ix, iy) < angiogenesis.threshold {
        return false;
    }
    let (gx, gy) = gradient(vegf, ix, iy);
    let norm = gx.hypot(gy);
    if !(norm.is_finite() && norm > 0.0) {
        return false;
    }
    let x = tip.x + angiogenesis.speed * gx / norm;
    let y = tip.y + angiogenesis.speed * gy / norm;
    // A disc laid at the edge would be cut off by it.
    let margin = f64::from(VESSEL_RADIUS);
    if x < margin || y < margin || x > f64::from(w - 1) - margin || y > f64::from(h - 1) - margin {
        return false;
    }
    tip.x = x;
    tip.y = y;
    let centre = (x.round() as u32, y.round() as u32);
    if tip.trail.last() == Some(&centre) {
        return true;
    }
    let clear = distance(centre, tip.origin) > CLEAR_OF_ORIGIN;
    let met = if clear {
        meets_another_vessel(lattice, centre, &tip.trail)
    } else {
        None
    };
    tree.lay(lattice, centre, tip.sprout);
    *laid = true;
    if let Some(pixel) = met {
        // A loop closes through both sprouts and whatever each grew from.
        tree.perfuse(tip.sprout);
        if let Some(other) = tree.owner(pixel) {
            tree.perfuse(other);
        }
    }
    tip.trail.push(centre);
    if tip.trail.len() > TRAIL_KEPT {
        tip.trail.remove(0);
    }
    met.is_none()
}

fn distance((ax, ay): (u32, u32), (bx, by): (u32, u32)) -> f64 {
    (f64::from(ax) - f64::from(bx)).hypot(f64::from(ay) - f64::from(by))
}

/// The first pixel of vessel, in reading order, that a disc laid at `centre`
/// would come within a pixel of and that is not one of the discs in `trail`.
fn meets_another_vessel(
    lattice: &CpmLattice,
    centre: (u32, u32),
    trail: &[(u32, u32)],
) -> Option<(u32, u32)> {
    let (w, h) = lattice.dims();
    let reach = VESSEL_RADIUS + 1;
    let radius = f64::from(VESSEL_RADIUS);
    let rows = centre.1.saturating_sub(reach)..(centre.1 + reach + 1).min(h);
    let cols = centre.0.saturating_sub(reach)..(centre.0 + reach + 1).min(w);
    rows.flat_map(|y| cols.clone().map(move |x| (x, y)))
        .find(|&p| {
            lattice.occupant(p.0, p.1) == OBSTACLE
                && distance(p, centre) <= f64::from(reach)
                && trail.iter().all(|&t| distance(p, t) > radius)
        })
}

/// The field's gradient at a pixel, by central differences, one sided at
/// the edge of the box.
fn gradient(field: &ScalarField, x: u32, y: u32) -> (f64, f64) {
    let (w, h) = field.dims();
    let (xl, xr) = (x.saturating_sub(1), (x + 1).min(w - 1));
    let (yl, yr) = (y.saturating_sub(1), (y + 1).min(h - 1));
    let dx = f64::from((xr - xl).max(1));
    let dy = f64::from((yr - yl).max(1));
    (
        (field.concentration_at(xr, y) - field.concentration_at(xl, y)) / dx,
        (field.concentration_at(x, yr) - field.concentration_at(x, yl)) / dy,
    )
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]

    use super::*;
    use cellweave_core::TimeStep;

    /// A box with a vessel at the centre and a field rising towards the
    /// right edge, as if the tissue there were calling.
    fn calling_from_the_right() -> (CpmLattice, ScalarField) {
        let mut lattice = CpmLattice::new(80, 40, 10.0).unwrap();
        lattice.add_obstacle_disc(Pos2 { x: 20, y: 20 }, 4);
        let mut field = ScalarField::new(80, 40, 0.2, 0.0005, 0.0).unwrap();
        for y in 0..40 {
            field.set_production(79, y, 0.5);
        }
        field.relax(TimeStep(1.0), 1e-9, 200_000).unwrap();
        (lattice, field)
    }

    fn vessel(lattice: &CpmLattice) -> usize {
        lattice.grid.iter().filter(|&&p| p == OBSTACLE).count()
    }

    #[test]
    fn a_vessel_sprouts_towards_the_call_and_lays_vessel_behind_its_tip() {
        let (mut lattice, field) = calling_from_the_right();
        let before = vessel(&lattice);
        let call = Angiogenesis {
            threshold: field.concentration_at(25, 20) * 0.9,
            speed: 0.5,
            reach: 20.0,
            perfusion: Perfusion::AtOnce,
        };
        let mut vessels = Vessels::default();
        for _ in 0..40 {
            vessels.grow(&mut lattice, &field, &call);
        }
        assert_eq!(
            vessels.sprouted(),
            1,
            "one wall reads the most, and spacing keeps the rest"
        );
        let tip = &vessels.tips()[0];
        assert!(tip.x > 35.0, "the tip climbed to the right: {}", tip.x);
        assert!(
            (tip.y - 20.0).abs() < 2.0,
            "and kept to the gradient: {}",
            tip.y
        );
        assert!(vessel(&lattice) > before + 40, "with vessel laid behind it");
    }

    #[test]
    fn nothing_sprouts_below_the_threshold() {
        let (mut lattice, field) = calling_from_the_right();
        let before = vessel(&lattice);
        let call = Angiogenesis {
            threshold: field.concentration_at(79, 20) * 2.0,
            speed: 0.5,
            reach: 20.0,
            perfusion: Perfusion::AtOnce,
        };
        let mut vessels = Vessels::default();
        for _ in 0..40 {
            assert!(!vessels.grow(&mut lattice, &field, &call));
        }
        assert_eq!(vessels.sprouted(), 0);
        assert_eq!(vessel(&lattice), before);
    }

    #[test]
    fn a_tip_that_reaches_another_vessel_joins_it_stops_and_carries_blood() {
        let (mut lattice, field) = calling_from_the_right();
        lattice.add_obstacle_disc(Pos2 { x: 40, y: 20 }, 3);
        let call = Angiogenesis {
            threshold: field.concentration_at(25, 20) * 0.9,
            speed: 0.5,
            reach: 20.0,
            perfusion: Perfusion::OnJoining,
        };
        let mut tree = Tree {
            sprouts: vec![Sprout {
                parent: None,
                perfused: false,
            }],
            owner: vec![None; 80 * 40],
            width: 80,
        };
        let mut tip = Tip::new((24, 20), 0);
        let mut laid = false;
        let mut steps = 0;
        while advance(&mut tip, &mut lattice, &field, &call, &mut tree, &mut laid) {
            steps += 1;
            assert!(steps < 100, "the tip never stopped");
            if lattice.occupant(30, 20) == OBSTACLE {
                assert!(!tree.carries_blood((30, 20)), "a blind end carries nothing");
            }
        }
        assert!(
            tip.x < 37.5,
            "it stopped at the other vessel's wall: {}",
            tip.x
        );
        assert_eq!(lattice.occupant(33, 20), OBSTACLE, "the two are joined");
        assert!(
            tree.carries_blood((30, 20)),
            "and blood flows through the loop"
        );
    }

    #[test]
    fn a_sprout_that_joins_nothing_carries_no_blood_unless_perfused_at_once() {
        for (perfusion, blood) in [(Perfusion::OnJoining, false), (Perfusion::AtOnce, true)] {
            let (mut lattice, field) = calling_from_the_right();
            let call = Angiogenesis {
                threshold: field.concentration_at(25, 20) * 0.9,
                speed: 0.5,
                reach: 20.0,
                perfusion,
            };
            let mut vessels = Vessels::default();
            for _ in 0..40 {
                vessels.grow(&mut lattice, &field, &call);
            }
            assert_eq!(
                lattice.occupant(32, 20),
                OBSTACLE,
                "the sprout reached here"
            );
            assert_eq!(vessels.carries_blood(32, 20), blood, "{perfusion:?}");
            assert!(
                vessels.carries_blood(20, 20),
                "the vessel laid at the start always does"
            );
        }
    }
}
