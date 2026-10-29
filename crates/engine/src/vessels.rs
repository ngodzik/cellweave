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
//! With a capillary lying in the plane, blood may instead flow as pressure
//! drives it: in at the capillary's left end, out at its right, through
//! every loop with a difference of pressure across it, by Poiseuille's law
//! (see `flow`). A blind end then carries nothing, and so does a loop closed
//! on a blind end.
//!
//! Nothing here is drawn at random. Which wall sprouts is the one reading the
//! most VEGF, and where a tip goes is where the gradient points, so a run
//! with vessels that may grow is the run without them, to the last bit, until
//! the first sprout lays its first pixel.

use std::collections::BTreeMap;

use cellweave_core::{Angiogenesis, DiffusionField, InTheWay, Perfusion, Pos2};

use crate::cpm::{CpmLattice, OBSTACLE};
use crate::flow::{self, Segment};
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

/// Flow through a vessel, as a share of what the capillary alone carries,
/// below which it counts as carrying no blood. Far below any real flow, far
/// above the rounding of the solve.
const FLOWING: f64 = 1e-6;

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

/// Where a sprout's tip joined another vessel.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Join {
    /// The sprout it joined, or `None` for the vessel the scenario laid.
    host: Option<usize>,
    /// Where along it: the index of the host's disc it met, or for the
    /// vessel the scenario laid, the column of the pixel it met.
    at: u32,
}

/// A capillary lying across the box, which blood enters at its left end and
/// leaves at its right.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) struct Capillary {
    /// The row it runs along.
    pub row: u32,
    /// Its length, the width of the box.
    pub width: u32,
    /// Half its width.
    pub radius: u32,
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
    /// Centre of every disc of vessel it laid.
    centres: Vec<(u32, u32)>,
    /// Whether it has regressed and its vessel been taken back.
    regressed: bool,
    /// Where its tip joined another vessel, if it has.
    joined: Option<Join>,
    /// The step since which it has carried no blood: when it started, or
    /// when the flow through it last stopped.
    dry_since: u64,
    /// When blood flows as pressure drives it, the flow through each of its
    /// discs, as a share of what the capillary alone carries.
    flow: Vec<f64>,
}

impl Sprout {
    fn new(parent: Option<usize>, perfused: bool, origin: (u32, u32), started: u64) -> Self {
        Self {
            parent,
            perfused,
            origin,
            centres: Vec::new(),
            regressed: false,
            joined: None,
            dry_since: started,
            flow: Vec::new(),
        }
    }

    /// Distance along it from where it started to each point of it: its
    /// origin, then the centre of each disc.
    fn arc(&self) -> Vec<f64> {
        let mut along = Vec::with_capacity(self.centres.len() + 1);
        let mut total = 0.0;
        let mut last = self.origin;
        along.push(0.0);
        for &c in &self.centres {
            total += distance(last, c);
            along.push(total);
            last = c;
        }
        along
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
    /// By pixel, for one a sprout laid: the index of the disc that laid it.
    disc: Vec<u32>,
    width: u32,
    /// The capillary blood flows through, when it flows as pressure drives
    /// it rather than by the rule of `Perfusion`.
    capillary: Option<Capillary>,
    /// Whether blood flows as pressure drives it.
    flowing: bool,
    /// Joins and regressions so far: the network blood flows through changes
    /// only with them.
    changes: u64,
    /// The value of `changes` the flow was last solved for.
    solved: u64,
    /// When blood flows as pressure drives it, the flow through the
    /// capillary at each column, as a share of what it alone carries.
    capillary_flow: Vec<f64>,
}

impl Tree {
    fn index(&self, (x, y): (u32, u32)) -> usize {
        y as usize * self.width as usize + x as usize
    }

    fn owner(&self, pixel: (u32, u32)) -> Option<usize> {
        self.owner.get(self.index(pixel)).copied().flatten()
    }

    fn disc_at(&self, pixel: (u32, u32)) -> u32 {
        self.disc.get(self.index(pixel)).copied().unwrap_or(0)
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
        let disc = self.sprouts[sprout].centres.len() as u32;
        for pixel in fresh {
            let i = self.index(pixel);
            self.owner[i] = Some(sprout);
            self.disc[i] = disc;
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
        for (i, (owner, disc)) in self.owner.iter_mut().zip(&mut self.disc).enumerate() {
            if !owner.is_some_and(|s| going[s]) {
                continue;
            }
            let pixel = ((i % width) as u32, (i / width) as u32);
            let keeper = self.sprouts.iter().enumerate().find_map(|(s, other)| {
                if other.regressed {
                    return None;
                }
                other
                    .centres
                    .iter()
                    .position(|&c| squared_distance(c, pixel) <= r2)
                    .map(|k| (s, k as u32))
            });
            *owner = keeper.map(|(s, _)| s);
            *disc = keeper.map_or(0, |(_, k)| k);
            if keeper.is_none() {
                lattice.clear_obstacle(pixel.0, pixel.1);
            }
        }
        self.changes += 1;
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
        match self.owner(pixel) {
            None => true,
            Some(s) if self.flowing => self.flow_at(s, self.disc_at(pixel)) > FLOWING,
            Some(s) => self.sprouts[s].perfused,
        }
    }

    /// The flow through the vessel at `pixel`, as a share of what the
    /// capillary alone carries: through the capillary at its column, or
    /// through the disc of a sprout that laid it. Zero unless blood flows as
    /// pressure drives it.
    fn flow_through(&self, pixel: (u32, u32)) -> f64 {
        match self.owner(pixel) {
            None => self
                .capillary_flow
                .get(pixel.0 as usize)
                .copied()
                .unwrap_or(0.0),
            Some(s) => self.flow_at(s, self.disc_at(pixel)),
        }
    }

    /// The flow through disc `disc` of `sprout`, as a share of what the
    /// capillary alone carries; zero unless blood flows as pressure drives it.
    fn flow_at(&self, sprout: usize, disc: u32) -> f64 {
        self.sprouts[sprout]
            .flow
            .get(disc as usize)
            .copied()
            .unwrap_or(0.0)
    }

    /// Solve the flow through the network, if it changed since the last
    /// solve, and with it which sprouts carry blood. `mcs` is the step being
    /// taken: a sprout the flow has just left is dry from then.
    fn solve_flow(&mut self, mcs: u64) {
        let Some(capillary) = self.capillary else {
            return;
        };
        if !self.flowing || self.changes == self.solved {
            return;
        }
        self.solved = self.changes;
        let (segments, placed, held, nodes) = self.network(capillary);
        let pressure = flow::pressures(nodes, &segments, &held);
        let q = flow::flows(&segments, &pressure);
        let alone = flow::conductance(
            f64::from(capillary.radius),
            f64::from(capillary.width.saturating_sub(1)),
        );
        for sprout in &mut self.sprouts {
            sprout.flow = vec![0.0; sprout.centres.len()];
        }
        self.capillary_flow = vec![0.0; capillary.width as usize];
        for (&(line, u, v), q) in placed.iter().zip(&q) {
            if line == 0 {
                let share = q.abs() / alone;
                for x in u.min(v)..=u.max(v) {
                    if let Some(f) = self.capillary_flow.get_mut(x as usize) {
                        *f = share;
                    }
                }
                continue;
            }
            // Position k along a sprout is its origin for 0 and its disc k - 1
            // after; a disc at a junction belongs to both segments it joins.
            let share = q.abs() / alone;
            let flow = &mut self.sprouts[line - 1].flow;
            for k in u.min(v).max(1)..=u.max(v) {
                if let Some(f) = flow.get_mut(k as usize - 1) {
                    *f = f.max(share);
                }
            }
        }
        for sprout in &mut self.sprouts {
            let wet = !sprout.regressed && sprout.flow.iter().any(|&f| f > FLOWING);
            if sprout.perfused && !wet {
                sprout.dry_since = mcs;
            }
            sprout.perfused = wet;
        }
    }

    /// The network: segments between junctions; for each segment, the line
    /// it lies along and between which positions; the junctions held at a
    /// pressure; and how many junctions there are. Line 0 is the capillary,
    /// where a position is a column, and line `s + 1` is sprout `s`, where
    /// position 0 is its origin and position `k + 1` its disc `k`.
    #[allow(clippy::type_complexity)]
    fn network(
        &self,
        capillary: Capillary,
    ) -> (
        Vec<Segment>,
        Vec<(usize, u32, u32)>,
        Vec<(usize, f64)>,
        usize,
    ) {
        let lines = self.sprouts.len() + 1;
        let mut at: Vec<BTreeMap<u32, usize>> = vec![BTreeMap::new(); lines];
        let mut nodes = 0;
        let mut same: Vec<(usize, usize)> = Vec::new();
        let mut node = |at: &mut Vec<BTreeMap<u32, usize>>, line: usize, pos: u32| -> usize {
            *at[line].entry(pos).or_insert_with(|| {
                nodes += 1;
                nodes - 1
            })
        };
        let line_of = |host: Option<usize>| host.map_or(0, |s| s + 1);
        let inlet = node(&mut at, 0, 0);
        let outlet = node(&mut at, 0, capillary.width.saturating_sub(1));
        let standing = |host: Option<usize>| host.is_none_or(|h| !self.sprouts[h].regressed);
        for (s, sprout) in self.sprouts.iter().enumerate() {
            if sprout.regressed || !standing(sprout.parent) {
                continue;
            }
            let start = match sprout.parent {
                None => sprout.origin.0,
                Some(_) => self.disc_at(sprout.origin) + 1,
            };
            let joint = node(&mut at, line_of(sprout.parent), start);
            if let Some(&mine) = at[s + 1].get(&0) {
                same.push((mine, joint));
            }
            at[s + 1].insert(0, joint);
            let end = sprout.centres.len() as u32;
            match sprout.joined {
                Some(j) if standing(j.host) && end > 0 => {
                    // On the capillary a join is at a column; on a sprout, at
                    // a disc, which is one position along from its origin.
                    let along = match j.host {
                        None => j.at,
                        Some(_) => j.at + 1,
                    };
                    let joint = node(&mut at, line_of(j.host), along);
                    if let Some(&mine) = at[s + 1].get(&end) {
                        same.push((mine, joint));
                    }
                    at[s + 1].insert(end, joint);
                }
                _ => {
                    node(&mut at, s + 1, end);
                }
            }
        }
        // Junctions that are one place, reached two ways, become one.
        let mut root: Vec<usize> = (0..nodes).collect();
        fn find(root: &mut [usize], mut n: usize) -> usize {
            while root[n] != n {
                root[n] = root[root[n]];
                n = root[n];
            }
            n
        }
        for (a, b) in same {
            let (ra, rb) = (find(&mut root, a), find(&mut root, b));
            if ra != rb {
                root[ra.max(rb)] = ra.min(rb);
            }
        }
        let mut segments = Vec::new();
        let mut placed = Vec::new();
        let capillary_radius = f64::from(capillary.radius);
        let sprout_radius = f64::from(VESSEL_RADIUS);
        for (line, positions) in at.iter().enumerate() {
            let arc = (line > 0).then(|| self.sprouts[line - 1].arc());
            let mut previous: Option<(u32, usize)> = None;
            for (&pos, &n) in positions {
                if let Some((p, m)) = previous {
                    let (length, radius) = match &arc {
                        None => (f64::from(pos - p), capillary_radius),
                        Some(arc) => (arc[pos as usize] - arc[p as usize], sprout_radius),
                    };
                    segments.push(Segment {
                        a: find(&mut root, m),
                        b: find(&mut root, n),
                        conductance: flow::conductance(radius, length),
                    });
                    placed.push((line, p, pos));
                }
                previous = Some((pos, n));
            }
        }
        let held = vec![
            (find(&mut root, inlet), 1.0),
            (find(&mut root, outlet), 0.0),
        ];
        (segments, placed, held, nodes)
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

    /// Let blood flow through `capillary`, in at its left end and out at
    /// its right, when the scenario says it flows as pressure drives it.
    pub(crate) fn lie_along(&mut self, capillary: Capillary, flowing: bool) {
        self.tree.capillary = Some(capillary);
        // Known from the start, not only from the first growth, since what
        // the capillary delivers is read before anything grows.
        self.tree.flowing = flowing;
        // Before any sprout, the capillary carries all of its own flow.
        self.tree.capillary_flow = vec![1.0; capillary.width as usize];
    }

    /// The flow through the vessel at `(x, y)`, as a share of what the
    /// capillary alone carries. Zero unless blood flows as pressure drives
    /// it, and meaningless for a pixel that is not vessel.
    pub(crate) fn flow_through(&self, x: u32, y: u32) -> f64 {
        if !self.tree.flowing {
            return 0.0;
        }
        self.tree.flow_through((x, y))
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
            self.tree.disc = vec![0; w as usize * h as usize];
            self.tree.width = w;
            self.tree.flowing = angiogenesis.perfusion == Perfusion::Flow;
        }
        let mut laid = false;
        let tree = &mut self.tree;
        self.tips
            .retain_mut(|tip| advance(tip, lattice, vegf, angiogenesis, tree, &mut laid));
        self.tree.solve_flow(mcs);
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

    /// Take back every sprout that has carried no blood for `after` steps or
    /// more by `mcs`, since it started or since the flow through it stopped,
    /// and with it every sprout grown from it, which carries none either:
    /// blood reaching one would have reached what it grew from.
    fn regress_dry(&mut self, lattice: &mut CpmLattice, after: u64, mcs: u64) {
        for i in 0..self.tree.sprouts.len() {
            let sprout = &self.tree.sprouts[i];
            let due = sprout.dry_since.saturating_add(after) <= mcs;
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
        let host = tree.owner(pixel);
        let at = match host {
            None => pixel.0,
            Some(_) => tree.disc_at(pixel),
        };
        tree.sprouts[tip.sprout].joined = Some(Join { host, at });
        tree.changes += 1;
        // A loop closes through both sprouts and whatever each grew from,
        // unless blood goes only where pressure drives it.
        if !tree.flowing {
            tree.perfuse(tip.sprout);
            if let Some(other) = host {
                tree.perfuse(other);
            }
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
            half_delivery_flow: None,
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
            half_delivery_flow: None,
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
            half_delivery_flow: None,
            in_the_way: InTheWay::Crushed,
        };
        let mut tree = Tree {
            sprouts: vec![Sprout::new(None, false, (24, 20), 0)],
            owner: vec![None; 80 * 40],
            disc: vec![0; 80 * 40],
            width: 80,
            ..Tree::default()
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
                half_delivery_flow: None,
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
            half_delivery_flow: None,
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
            disc: vec![0; 80 * 40],
            width: 80,
            ..Tree::default()
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
                half_delivery_flow: None,
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

    /// A capillary along row 10 of an 80 by 40 box, a bypass from column 20
    /// back to column 40, a blind sprout down from column 60, and a loop from
    /// the blind sprout back onto itself.
    fn a_network_laid_by_hand() -> (CpmLattice, Tree) {
        let mut lattice = CpmLattice::new(80, 40, 10.0).unwrap();
        for x in 0..80 {
            lattice.add_obstacle_disc(Pos2 { x, y: 10 }, 2);
        }
        let mut tree = Tree {
            sprouts: Vec::new(),
            owner: vec![None; 80 * 40],
            disc: vec![0; 80 * 40],
            width: 80,
            capillary: Some(Capillary {
                row: 10,
                width: 80,
                radius: 2,
            }),
            flowing: true,
            changes: 0,
            solved: 0,
            capillary_flow: vec![1.0; 80],
        };
        let sprout =
            |tree: &mut Tree, lattice: &mut CpmLattice, origin: (u32, u32), path: &[(u32, u32)]| {
                let parent = tree.owner(origin);
                tree.sprouts.push(Sprout::new(parent, false, origin, 0));
                let s = tree.sprouts.len() - 1;
                for &c in path {
                    tree.lay(lattice, c, s, None);
                }
                s
            };
        let line = |from: (u32, u32), to: (u32, u32)| -> Vec<(u32, u32)> {
            let n = from.0.abs_diff(to.0).max(from.1.abs_diff(to.1));
            (1..=n)
                .map(|k| {
                    let t = f64::from(k) / f64::from(n);
                    let x = f64::from(from.0) + t * (f64::from(to.0) - f64::from(from.0));
                    let y = f64::from(from.1) + t * (f64::from(to.1) - f64::from(from.1));
                    (x.round() as u32, y.round() as u32)
                })
                .collect()
        };
        let bypass: Vec<(u32, u32)> = [
            line((20, 12), (20, 25)),
            line((20, 25), (40, 25)),
            line((40, 25), (40, 15)),
        ]
        .concat();
        let s = sprout(&mut tree, &mut lattice, (20, 12), &bypass);
        tree.sprouts[s].joined = Some(Join { host: None, at: 40 });
        let blind = sprout(&mut tree, &mut lattice, (60, 12), &line((60, 12), (60, 30)));
        let loop_back: Vec<(u32, u32)> = [
            line((62, 18), (70, 18)),
            line((70, 18), (70, 26)),
            line((70, 26), (65, 26)),
        ]
        .concat();
        let s = sprout(&mut tree, &mut lattice, (62, 18), &loop_back);
        let at = tree.disc_at((62, 26));
        tree.sprouts[s].joined = Some(Join {
            host: Some(blind),
            at,
        });
        tree.changes += 1;
        (lattice, tree)
    }

    #[test]
    fn blood_flows_through_a_bypass_and_not_through_a_blind_end() {
        let (_, mut tree) = a_network_laid_by_hand();
        tree.solve_flow(0);
        assert!(tree.carries_blood((5, 10)), "the capillary itself");
        assert!(tree.carries_blood((30, 25)), "the bypass");
        assert!(!tree.carries_blood((60, 28)), "the blind sprout");
        assert!(
            !tree.carries_blood((70, 22)),
            "a loop closed on a blind end has no difference of pressure across it"
        );
        let through = tree.flow_at(0, 20);
        assert!(
            through > 0.0 && through < 1.0,
            "the bypass takes a share of the capillary's flow: {through}"
        );
        let perfused: Vec<bool> = tree.sprouts.iter().map(|s| s.perfused).collect();
        assert_eq!(perfused, [true, false, false]);
    }

    #[test]
    fn a_narrower_vessel_takes_less_of_the_flow() {
        // The same bypass, against a capillary twice as wide: its share of
        // the flow shrinks, by Poiseuille's fourth power.
        let (_, mut narrow) = a_network_laid_by_hand();
        narrow.solve_flow(0);
        let (_, mut wide) = a_network_laid_by_hand();
        wide.capillary = wide.capillary.map(|c| Capillary { radius: 4, ..c });
        wide.solve_flow(0);
        assert!(wide.flow_at(0, 20) < narrow.flow_at(0, 20) / 4.0);
    }

    #[test]
    fn a_sprout_the_flow_leaves_is_dry_from_then() {
        let (mut lattice, mut tree) = a_network_laid_by_hand();
        tree.solve_flow(5);
        assert!(tree.sprouts[0].perfused);
        // Taking back the capillary end it joined at: without it the bypass
        // is a blind end, and stops carrying blood at step 9.
        tree.sprouts[0].joined = None;
        tree.changes += 1;
        tree.solve_flow(9);
        assert!(!tree.sprouts[0].perfused);
        assert_eq!(tree.sprouts[0].dry_since, 9);
        let mut vessels = Vessels {
            tips: Vec::new(),
            tree,
        };
        vessels.regress_dry(&mut lattice, 10, 18);
        assert!(!vessels.tree.sprouts[0].regressed, "dry for 9 steps only");
        vessels.regress_dry(&mut lattice, 10, 19);
        assert!(vessels.tree.sprouts[0].regressed);
    }

    #[test]
    fn a_join_at_a_blind_sprouts_first_disc_makes_it_carry_blood() {
        // A blind sprout down from column 60, and a tip from column 75 that
        // comes back along to meet its first disc: blood runs from column
        // 60 down that first disc, through the second sprout, to column 75.
        let (mut lattice, mut tree) = a_network_laid_by_hand();
        tree.sprouts.truncate(2);
        let blind = 1;
        tree.sprouts.push(Sprout::new(None, false, (75, 12), 0));
        let s = tree.sprouts.len() - 1;
        for c in [
            (75, 14),
            (75, 15),
            (74, 15),
            (72, 15),
            (70, 15),
            (68, 15),
            (66, 15),
            (64, 15),
        ] {
            tree.lay(&mut lattice, c, s, None);
        }
        let first = tree.sprouts[blind].centres[0];
        let at = tree.disc_at(first);
        assert_eq!(
            at, 0,
            "the pixel met was laid by the blind sprout's first disc"
        );
        tree.sprouts[s].joined = Some(Join {
            host: Some(blind),
            at,
        });
        tree.changes += 1;
        tree.solve_flow(0);
        assert!(
            tree.flow_at(blind, 0) > FLOWING,
            "its first disc carries blood"
        );
        assert!(tree.sprouts[blind].perfused);
        assert!(tree.flow_at(blind, 10) <= FLOWING, "and beyond it, nothing");
    }

    #[test]
    fn the_capillary_carries_less_where_a_bypass_takes_its_share() {
        let (_, mut tree) = a_network_laid_by_hand();
        tree.capillary_flow = vec![1.0; 80];
        tree.solve_flow(0);
        let (before, beside, after) = (
            tree.flow_through((10, 10)),
            tree.flow_through((30, 10)),
            tree.flow_through((50, 10)),
        );
        assert!((before - after).abs() < 1e-9, "what goes in comes out");
        assert!(before > 1.0, "a bypass lowers the resistance: {before}");
        let bypass = tree.flow_through((30, 25));
        assert!(
            (beside + bypass - before).abs() < 1e-9,
            "the bypassed stretch and the bypass share it: {beside} + {bypass} against {before}"
        );
    }
}
