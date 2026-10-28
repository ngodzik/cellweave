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
//! grew from. When the scenario says so, a sprout that has carried no blood
//! for long enough regresses, with every sprout grown from it: its vessel is
//! taken back, and the wall it started from may sprout again.
//!
//! Nothing here is drawn at random. Which wall sprouts is the one reading the
//! most VEGF, and where a tip goes is where the gradient points, so a run
//! with vessels that may grow is the run without them, to the last bit, until
//! the first sprout lays its first pixel.

use cellweave_core::{Angiogenesis, DiffusionField, InTheWay, Perfusion, Pos2};

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
    /// The wall pixel it started from.
    origin: (u32, u32),
    /// The step it started at.
    started: u64,
    /// Centre of every disc of vessel it laid.
    centres: Vec<(u32, u32)>,
    /// Whether it has regressed and its vessel been taken back.
    regressed: bool,
}

impl Sprout {
    fn new(parent: Option<usize>, perfused: bool, origin: (u32, u32), started: u64) -> Self {
        Self {
            parent,
            perfused,
            origin,
            started,
            centres: Vec::new(),
            regressed: false,
        }
    }
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
    /// pixels it made vessel. With a heading, the way the tip is moving, the
    /// cells in the disc are pushed aside; without one, crushed.
    fn lay(
        &mut self,
        lattice: &mut CpmLattice,
        centre: (u32, u32),
        sprout: usize,
        push: Option<(f64, f64)>,
    ) {
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
        let at = Pos2 {
            x: centre.0,
            y: centre.1,
        };
        match push {
            Some(heading) => {
                lattice.push_obstacle_disc(at, r, heading);
            }
            None => lattice.add_obstacle_disc(at, r),
        }
        for pixel in fresh {
            let i = self.index(pixel);
            self.owner[i] = Some(sprout);
        }
        self.sprouts[sprout].centres.push(centre);
    }

    /// Take back the vessel of `root` and of every sprout grown from it.
    /// A pixel one of them laid is given back to the tissue, unless it lies
    /// in a disc a sprout still standing laid too, which then keeps it.
    fn regress(&mut self, lattice: &mut CpmLattice, root: usize) {
        // A sprout comes after the one it grew from, so one pass in order
        // finds every sprout grown from the root, however indirectly.
        let mut going = vec![false; self.sprouts.len()];
        for (i, sprout) in self.sprouts.iter().enumerate().skip(root) {
            going[i] = i == root || sprout.parent.is_some_and(|p| going[p]);
        }
        for (sprout, &gone) in self.sprouts.iter_mut().zip(&going) {
            sprout.regressed |= gone;
        }
        let r2 = i64::from(VESSEL_RADIUS * VESSEL_RADIUS);
        let width = self.width as usize;
        for (i, owner) in self.owner.iter_mut().enumerate() {
            if !owner.is_some_and(|s| going[s]) {
                continue;
            }
            let pixel = ((i % width) as u32, (i / width) as u32);
            let keeper = self.sprouts.iter().position(|other| {
                !other.regressed
                    && other
                        .centres
                        .iter()
                        .any(|&c| squared_distance(c, pixel) <= r2)
            });
            *owner = keeper;
            if keeper.is_none() {
                lattice.clear_obstacle(pixel.0, pixel.1);
            }
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

/// Every sprout still growing, every sprout ever started, and which of them
/// carry blood.
#[derive(Debug, Clone, Default)]
pub(crate) struct Vessels {
    tips: Vec<Tip>,
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

    /// How many sprouts have started, growing or not, standing or not.
    pub(crate) fn sprouted(&self) -> usize {
        self.tree.sprouts.len()
    }

    /// How many sprouts have regressed.
    pub(crate) fn regressed(&self) -> usize {
        self.tree.sprouts.iter().filter(|s| s.regressed).count()
    }

    /// Move every tip, then take back every sprout that has gone too long
    /// without blood, then let the wall reading the most VEGF sprout, at most
    /// once. `mcs` is the step being taken. Returns whether any vessel was
    /// laid, which is when cells may have lost pixels to it.
    pub(crate) fn grow(
        &mut self,
        lattice: &mut CpmLattice,
        vegf: &ScalarField,
        angiogenesis: &Angiogenesis,
        mcs: u64,
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
        // Checked after the tips move, so that a sprout joining another
        // vessel on its last step stands.
        if let Some(after) = angiogenesis.regress_after {
            self.regress_dry(lattice, after, mcs);
        }
        if let Some(wall) = self.strongest_free_wall(lattice, vegf, angiogenesis.threshold) {
            let sprout = self.tree.sprouts.len();
            self.tree.sprouts.push(Sprout::new(
                self.tree.owner(wall),
                angiogenesis.perfusion == Perfusion::AtOnce,
                wall,
                mcs,
            ));
            self.tips.push(Tip::new(wall, sprout));
        }
        laid
    }

    /// Take back every sprout that started `after` steps or more before
    /// `mcs` and carries no blood, and with it every sprout grown from it,
    /// which carries none either: blood reaching one would have reached
    /// what it grew from.
    fn regress_dry(&mut self, lattice: &mut CpmLattice, after: u64, mcs: u64) {
        for i in 0..self.tree.sprouts.len() {
            let sprout = &self.tree.sprouts[i];
            let due = sprout.started.saturating_add(after) <= mcs;
            if due && !sprout.perfused && !sprout.regressed {
                self.tree.regress(lattice, i);
            }
        }
        let sprouts = &self.tree.sprouts;
        self.tips.retain(|tip| !sprouts[tip.sprout].regressed);
    }

    /// The vessel pixel open to the tissue that reads the most VEGF, at
    /// least `threshold`, and is far enough from every tip and from where
    /// every sprout still standing started.
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
                    || self.tree.sprouts.iter().any(|s| {
                        !s.regressed && near(x, y, (f64::from(s.origin.0), f64::from(s.origin.1)))
                    });
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
    let push = (angiogenesis.in_the_way == InTheWay::Pushed).then(|| {
        let from = tip.trail.last().copied().unwrap_or(tip.origin);
        (
            f64::from(centre.0) - f64::from(from.0),
            f64::from(centre.1) - f64::from(from.1),
        )
    });
    tree.lay(lattice, centre, tip.sprout, push);
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

fn squared_distance((ax, ay): (u32, u32), (bx, by): (u32, u32)) -> i64 {
    let (dx, dy) = (i64::from(ax) - i64::from(bx), i64::from(ay) - i64::from(by));
    dx * dx + dy * dy
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
            regress_after: None,
            in_the_way: InTheWay::Crushed,
        };
        let mut vessels = Vessels::default();
        for mcs in 0..40 {
            vessels.grow(&mut lattice, &field, &call, mcs);
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
            regress_after: None,
            in_the_way: InTheWay::Crushed,
        };
        let mut vessels = Vessels::default();
        for mcs in 0..40 {
            assert!(!vessels.grow(&mut lattice, &field, &call, mcs));
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
            regress_after: None,
            in_the_way: InTheWay::Crushed,
        };
        let mut tree = Tree {
            sprouts: vec![Sprout::new(None, false, (24, 20), 0)],
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
                regress_after: None,
                in_the_way: InTheWay::Crushed,
            };
            let mut vessels = Vessels::default();
            for mcs in 0..40 {
                vessels.grow(&mut lattice, &field, &call, mcs);
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

    /// A sprout that never joins anything regresses `regress_after` steps
    /// after it started, and the wall it left may sprout again.
    #[test]
    fn a_blind_sprout_regresses_and_gives_its_pixels_back() {
        let (mut lattice, field) = calling_from_the_right();
        let before = lattice.grid.clone();
        let call = Angiogenesis {
            threshold: field.concentration_at(25, 20) * 0.9,
            speed: 0.5,
            reach: 20.0,
            perfusion: Perfusion::OnJoining,
            regress_after: Some(20),
            in_the_way: InTheWay::Crushed,
        };
        let mut vessels = Vessels::default();
        for mcs in 0..20 {
            vessels.grow(&mut lattice, &field, &call, mcs);
        }
        assert!(vessel(&lattice) > vessel_count(&before), "the sprout grew");
        assert_eq!(vessels.regressed(), 0, "and stands until its time is up");
        vessels.grow(&mut lattice, &field, &call, 20);
        assert_eq!(vessels.regressed(), 1);
        assert_eq!(
            lattice.grid, before,
            "every pixel it laid is given back, and only those"
        );
        assert_eq!(vessels.sprouted(), 2, "the wall it left sprouted again");
        assert_eq!(vessels.tips().len(), 1, "and only the new tip grows");
    }

    fn vessel_count(grid: &ndarray::Array2<u32>) -> usize {
        grid.iter().filter(|&&p| p == OBSTACLE).count()
    }

    /// A tree laid by hand: sprout 0 to the right of the vessel and sprout
    /// 1 grown from it, neither carrying blood; sprout 2 above the vessel and
    /// sprout 3 below sprout 0, overlapping it, both carrying blood.
    fn a_tree_laid_by_hand() -> (CpmLattice, Vessels) {
        let mut lattice = CpmLattice::new(80, 40, 10.0).unwrap();
        lattice.add_obstacle_disc(Pos2 { x: 20, y: 20 }, 4);
        let mut tree = Tree {
            sprouts: vec![
                Sprout::new(None, false, (24, 20), 0),
                Sprout::new(Some(0), false, (29, 18), 15),
                Sprout::new(None, true, (20, 16), 0),
                Sprout::new(None, true, (24, 22), 0),
            ],
            owner: vec![None; 80 * 40],
            width: 80,
        };
        let discs: [(usize, &[(u32, u32)]); 4] = [
            (0, &[(26, 20), (27, 20), (28, 20), (29, 20), (30, 20)]),
            (1, &[(29, 16), (29, 15), (29, 14)]),
            (2, &[(20, 13), (20, 12), (20, 11)]),
            (3, &[(26, 24), (27, 24)]),
        ];
        for (sprout, centres) in discs {
            for &centre in centres {
                tree.lay(&mut lattice, centre, sprout, None);
            }
        }
        (
            lattice,
            Vessels {
                tips: Vec::new(),
                tree,
            },
        )
    }

    #[test]
    fn a_sprout_regresses_with_what_grew_from_it_and_leaves_what_another_laid() {
        let (mut lattice, mut vessels) = a_tree_laid_by_hand();
        let laid = lattice.grid.clone();
        vessels.regress_dry(&mut lattice, 20, 19);
        assert_eq!(lattice.grid, laid, "nothing is due before step 20");

        vessels.regress_dry(&mut lattice, 20, 20);
        let regressed: Vec<bool> = vessels.tree.sprouts.iter().map(|s| s.regressed).collect();
        assert_eq!(
            regressed,
            [true, true, false, false],
            "sprout 1 goes with sprout 0 although its own time is not up"
        );
        assert_eq!(
            lattice.occupant(30, 20),
            0,
            "sprout 0's vessel is given back"
        );
        assert_eq!(lattice.occupant(29, 14), 0, "and sprout 1's");
        assert_eq!(lattice.occupant(20, 11), OBSTACLE, "sprout 2 carries blood");
        assert_eq!(
            lattice.occupant(27, 22),
            OBSTACLE,
            "sprout 0 laid this pixel first, inside a disc sprout 3 laid too"
        );
        assert_eq!(vessels.tree.owner((27, 22)), Some(3), "which now holds it");
        assert_eq!(
            lattice.occupant(23, 20),
            OBSTACLE,
            "the vessel the scenario laid never regresses"
        );
        assert!(vessels.carries_blood(27, 22));
    }

    /// A small cell sits in the path of the sprout: crushed, it is left with
    /// nothing; pushed, it keeps every pixel it had.
    #[test]
    fn a_cell_in_the_way_is_crushed_or_pushed_aside() {
        for (in_the_way, kept) in [(InTheWay::Crushed, false), (InTheWay::Pushed, true)] {
            let (mut lattice, field) = calling_from_the_right();
            let cell = lattice.add_cell(cellweave_core::CellKind::Tumor, Pos2 { x: 30, y: 20 }, 1);
            let size = lattice.volume(cell).unwrap().0;
            let call = Angiogenesis {
                threshold: field.concentration_at(25, 20) * 0.9,
                speed: 0.5,
                reach: 20.0,
                perfusion: Perfusion::AtOnce,
                regress_after: None,
                in_the_way,
            };
            let mut vessels = Vessels::default();
            for mcs in 0..40 {
                vessels.grow(&mut lattice, &field, &call, mcs);
            }
            assert_eq!(
                lattice.occupant(32, 20),
                OBSTACLE,
                "the sprout went through"
            );
            let left = lattice.volume(cell).unwrap().0;
            assert_eq!(
                left == size,
                kept,
                "{in_the_way:?}: {left} of {size} pixels"
            );
            if !kept {
                assert_eq!(left, 0, "crushed, it has nothing left");
            }
        }
    }
}
