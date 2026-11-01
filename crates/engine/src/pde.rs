//! Scalar extracellular diffusion field, explicit finite differences.
//!
//! Solves ∂c/∂t = D·∇²c − λ·c + S on a grid, flat or in a volume: a volume
//! stacks slices, and each voxel has two neighbours more. A flat field
//! visits and sums what it always did, in the same order.
//!
//! The walls of the simulated square let nothing through. Everything that arrives
//! comes from [`ScalarField::hold_bath`], and that is deliberate: when the edge of
//! the box is the only supply, the distance from the tissue to that edge decides
//! how much reaches it, and the size of the box becomes a biological parameter. A
//! bath held against the surface of the tissue removes the box from the answer,
//! and leaves the gradient to form where it belongs, inside the tissue, out of the
//! balance between what diffuses in and what the cells take up.
//!
//! Oxygen crosses a cell in a fraction of a second, a signalling network answers
//! in minutes, a cell cycle takes a day. At the scale where a cell decides
//! anything, the field is therefore already at equilibrium with the layout of the
//! cells. [`ScalarField::relax`] is how that is honoured: the field is iterated to
//! a steady state rather than advanced by some number of sub steps per lattice
//! step, which would leave it trailing the cells by an amount set by a knob with
//! no physical meaning.

use cellweave_core::{CellweaveError, DiffusionField, TimeStep};
use rayon::prelude::*;

/// Scalar diffusion field with no-flux walls, flat or in a volume.
///
/// Values are stored slice after slice, row after row, so that the first
/// slice is a flat field's every pixel and `(x, y)` alone names a pixel of it.
pub struct ScalarField {
    width: u32,
    height: u32,
    /// Slices: one for a flat field.
    depth: u32,
    /// Current concentration on each pixel.
    values: Vec<f64>,
    /// Scratch buffer for the next time step.
    next: Vec<f64>,
    /// Diffusion coefficient (pixels² per unit time).
    d: f64,
    /// First-order decay rate (per unit time).
    decay: f64,
    /// Uniform source term (concentration per unit time).
    source: f64,
    /// Pixels the bath holds, which `step` leaves alone.
    held: Vec<bool>,
    /// Sink per pixel, what the tissue takes up per unit time.
    uptake: Vec<f64>,
    /// Source per pixel, what the tissue secretes per unit time.
    production: Vec<f64>,
    /// Scratch: breadth-first frontier for the bath fill, by flat index.
    frontier: Vec<usize>,
}

/// A box of pixels, each side an inclusive range.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Span {
    x: (u32, u32),
    y: (u32, u32),
    z: (u32, u32),
}

/// Rows one parallel task takes at least.
///
/// A single row is too little work to hand to another thread, and too many
/// rows per task leave threads idle on a 200 pixel grid, where a half sweep
/// is only a few dozen tasks. Measured on the spheroid, four was the fastest of
/// one, four, eight, sixteen and thirty-two, with eight a close second.
const ROWS_PER_TASK: usize = 4;

/// Largest `D · dt` an explicit five point stencil stays stable at, on a unit
/// grid: one over twice the number of dimensions, so a third less in a volume.
const DIFFUSION_LIMIT: f64 = 0.25;

impl ScalarField {
    /// Create a flat field holding `initial` everywhere.
    ///
    /// # Errors
    ///
    /// [`CellweaveError::Config`] when either dimension is zero, or when the
    /// diffusion coefficient or the decay rate is negative.
    pub fn new(
        width: u32,
        height: u32,
        d: f64,
        decay: f64,
        initial: f64,
    ) -> Result<Self, CellweaveError> {
        Self::with_depth(width, height, 1, d, decay, initial)
    }

    /// Create a field `depth` slices deep holding `initial` everywhere: flat
    /// for a depth of one.
    ///
    /// # Errors
    ///
    /// As [`ScalarField::new`], and when the depth is zero.
    pub fn with_depth(
        width: u32,
        height: u32,
        depth: u32,
        d: f64,
        decay: f64,
        initial: f64,
    ) -> Result<Self, CellweaveError> {
        if width == 0 || height == 0 || depth == 0 {
            return Err(CellweaveError::Config(format!(
                "field dimensions must be positive, got {width}x{height}x{depth}"
            )));
        }
        if !d.is_finite() || !decay.is_finite() || d < 0.0 || decay < 0.0 {
            return Err(CellweaveError::Config(format!(
                "diffusion and decay must be finite and not negative, got D = {d}, decay = {decay}"
            )));
        }
        let n = width as usize * height as usize * depth as usize;
        Ok(Self {
            width,
            height,
            depth,
            values: vec![initial; n],
            next: vec![initial; n],
            d,
            decay,
            source: 0.0,
            held: vec![false; n],
            uptake: vec![0.0; n],
            production: vec![0.0; n],
            frontier: Vec::new(),
        })
    }

    /// Set a spatially uniform source term.
    #[must_use]
    pub fn with_source(mut self, source: f64) -> Self {
        self.source = source;
        self
    }

    /// Set how much one pixel takes up per unit time.
    ///
    /// This is how tissue consumes: as a sink in the equation, applied on every
    /// step, rather than an amount subtracted once, which relaxation would simply
    /// diffuse away. Coordinates outside the grid are ignored, since callers
    /// iterate over a lattice of the same size.
    pub fn set_uptake(&mut self, x: u32, y: u32, rate: f64) {
        if x >= self.width || y >= self.height {
            return;
        }
        let i = self.idx(x, y);
        self.uptake[i] = rate.max(0.0);
    }

    /// Set how much one pixel secretes per unit time.
    ///
    /// The mirror of [`ScalarField::set_uptake`]: a source in the equation,
    /// on every step, rather than an amount added once. Coordinates outside
    /// the grid are ignored.
    pub fn set_production(&mut self, x: u32, y: u32, rate: f64) {
        if x >= self.width || y >= self.height {
            return;
        }
        let i = self.idx(x, y);
        self.production[i] = rate.max(0.0);
    }

    /// Set how much the pixel at flat index `i` takes up per unit time, as
    /// [`ScalarField::set_uptake`] does by coordinates. An index past the end
    /// is ignored.
    pub fn set_uptake_at(&mut self, i: usize, rate: f64) {
        if let Some(u) = self.uptake.get_mut(i) {
            *u = rate.max(0.0);
        }
    }

    /// Set how much the pixel at flat index `i` secretes per unit time. An
    /// index past the end is ignored.
    pub fn set_production_at(&mut self, i: usize, rate: f64) {
        if let Some(p) = self.production.get_mut(i) {
            *p = rate.max(0.0);
        }
    }

    /// The concentration at flat index `i`, counted along rows, then down,
    /// then through the slices; zero past the end.
    pub fn value_at(&self, i: usize) -> f64 {
        self.values.get(i).copied().unwrap_or(0.0)
    }

    /// Pixels in the field, every slice counted.
    pub fn len(&self) -> usize {
        self.values.len()
    }

    /// Whether the field has no pixel, which construction never allows.
    pub fn is_empty(&self) -> bool {
        self.values.is_empty()
    }

    /// Width, height and depth.
    pub fn dims3(&self) -> (u32, u32, u32) {
        (self.width, self.height, self.depth)
    }

    /// Forget every source, before what the tissue secretes is read in again.
    pub fn clear_production(&mut self) {
        self.production.fill(0.0);
    }

    /// Set every pixel to `value`, as when what supplies a field is cut off
    /// and everything it held has been taken up.
    pub fn fill(&mut self, value: f64) {
        self.values.fill(value);
        self.next.fill(value);
    }

    /// Forget every sink, before the layout of the tissue is read in again.
    pub fn clear_uptake(&mut self) {
        self.uptake.fill(0.0);
    }

    /// Solve for the steady state directly, and say how many sweeps it took.
    ///
    /// This does not step forward in time. It sweeps the grid in red-black order,
    /// setting each pixel to the value that balances diffusion against decay,
    /// source and uptake given its neighbours as they stand, pushed a little
    /// past that value to hurry the slow modes along: successive over
    /// relaxation, with the result floored at zero since tissue cannot take up
    /// what is not there. Stepping in time would reach the same state, but the
    /// slowest mode of a region of size L decays in about L² steps, where a
    /// sweep like this needs about L.
    ///
    /// Converged means no pixel moved by more than `tolerance` over the last
    /// sweep. Pick it well below the uptake, since the per sweep change is what
    /// has not yet been balanced. `max_sweeps` is a hard limit that turns a
    /// field which will not settle into an error rather than a silent partial
    /// answer.
    ///
    /// The time step is unused and kept so callers that switch between stepping
    /// and relaxing keep the same shape; a steady state has no time step.
    ///
    /// # Errors
    ///
    /// [`CellweaveError::Simulation`] when `max_sweeps` pass without convergence.
    pub fn relax(
        &mut self,
        _dt: TimeStep,
        tolerance: f64,
        max_sweeps: usize,
    ) -> Result<usize, CellweaveError> {
        // A held pixel is a boundary condition and never changes, so only the
        // smallest box around the free ones is swept. A bath holds everything
        // outside the tissue, and a spheroid starts as a speck in a box sized
        // for the end of the run. Outside that box the two buffers have to
        // agree, since a sweep no longer writes there before they are swapped.
        let Some(span) = self.free_span() else {
            // Nothing to balance: the first sweep would move nothing.
            return Ok(1);
        };

        // The optimal factor for a Poisson problem on an N by N grid is close to
        // 2 / (1 + sin(π / N)). The plane stays a little under it, capped at
        // 1.95, which is what its runs were made with. A volume takes the
        // optimum for the box it sweeps: a factor that stays below two keeps
        // the floored iteration converging, and capped at 1.95 a box 250 on a
        // side settles about ten times slower.
        let omega = if self.depth > 1 {
            let longest = [span.x, span.y, span.z]
                .iter()
                .map(|&(lo, hi)| hi - lo + 1)
                .max()
                .unwrap_or(1);
            2.0 / (1.0 + (std::f64::consts::PI / f64::from(longest.max(2))).sin())
        } else {
            let n = f64::from(self.width.max(self.height));
            (2.0 / (1.0 + (std::f64::consts::PI / n).sin())).min(1.95)
        };
        self.next.copy_from_slice(&self.values);
        for done in 1..=max_sweeps {
            let moved = self.sor_sweep(omega, span);
            if moved < tolerance {
                return Ok(done);
            }
        }
        Err(CellweaveError::Simulation(format!(
            "diffusion did not settle within {max_sweeps} sweeps at tolerance {tolerance}"
        )))
    }

    /// One sweep in two halves, returning the largest change made to any pixel.
    ///
    /// The grid is coloured like a chessboard. Every neighbour of a red pixel is
    /// black, so all the red pixels can be balanced at once from the black ones
    /// as they stand, then all the black ones from the new red ones. That is the
    /// red-black ordering of the same over-relaxation, and it is what lets the
    /// grid be swept in parallel: no pixel ever reads a neighbour being written.
    fn sor_sweep(&mut self, omega: f64, span: Span) -> f64 {
        let red = self.half_sweep(omega, 0, span);
        let black = self.half_sweep(omega, 1, span);
        red.max(black)
    }

    /// The smallest box holding every pixel not held, or `None` when every
    /// pixel is.
    fn free_span(&self) -> Option<Span> {
        let (w, h) = (self.width as usize, self.height as usize);
        let mut span: Option<Span> = None;
        for (i, _) in self.held.iter().enumerate().filter(|(_, held)| !**held) {
            let (x, y, z) = ((i % w) as u32, ((i / w) % h) as u32, (i / (w * h)) as u32);
            let s = span.get_or_insert(Span {
                x: (x, x),
                y: (y, y),
                z: (z, z),
            });
            s.x = (s.x.0.min(x), s.x.1.max(x));
            s.y = (s.y.0.min(y), s.y.1.max(y));
            s.z.1 = z;
        }
        span
    }

    /// Balance every pixel of one colour, row by row in parallel.
    ///
    /// Reads `values` and writes `next`, then swaps them, so each pixel's new
    /// value is a function of the old grid alone. Which thread computes which
    /// row cannot change any value, and the largest change is a maximum, which
    /// does not depend on the order it is taken in: one thread or twenty, the
    /// field is the same to the last bit.
    fn half_sweep(&mut self, omega: f64, colour: u32, span: Span) -> f64 {
        let (width, height, depth) = (self.width, self.height, self.depth);
        let (d, decay, source) = (self.d, self.decay, self.source);
        let (values, held, uptake) = (&self.values, &self.held, &self.uptake);
        let production = &self.production;
        let (w, h) = (width as usize, height as usize);
        let slice = w * h;
        // Rows are counted down the stacked slices. From the span's first row
        // to its last, which in a volume takes in rows of the slices between
        // that lie outside the span and are skipped.
        let first = span.z.0 as usize * h + span.y.0 as usize;
        let last = span.z.1 as usize * h + span.y.1 as usize;

        let moved = self.next[first * w..(last + 1) * w]
            .par_chunks_mut(w)
            .with_min_len(ROWS_PER_TASK)
            .enumerate()
            .map(|(k, row)| {
                let stacked = first + k;
                let (y, z) = ((stacked % h) as u32, (stacked / h) as u32);
                let mut moved: f64 = 0.0;
                if y < span.y.0 || y > span.y.1 {
                    return moved;
                }
                for x in span.x.0..=span.x.1 {
                    let out = &mut row[x as usize];
                    let i = stacked * w + x as usize;
                    let before = values[i];
                    if (x + y + z) % 2 != colour || held[i] {
                        *out = before;
                        continue;
                    }
                    // Real neighbours only; a missing one at the wall is a
                    // mirror of this pixel, which the denominator accounts for.
                    let mut sum = 0.0;
                    let mut real = 0.0;
                    if x > 0 {
                        sum += values[i - 1];
                        real += 1.0;
                    }
                    if x + 1 < width {
                        sum += values[i + 1];
                        real += 1.0;
                    }
                    if y > 0 {
                        sum += values[i - w];
                        real += 1.0;
                    }
                    if y + 1 < height {
                        sum += values[i + w];
                        real += 1.0;
                    }
                    // Last, so that a flat field adds what it always did in
                    // the order it always did.
                    if z > 0 {
                        sum += values[i - slice];
                        real += 1.0;
                    }
                    if z + 1 < depth {
                        sum += values[i + slice];
                        real += 1.0;
                    }
                    // A pixel that secretes nothing adds an exact zero, so a
                    // field without production balances as it always did.
                    let balanced =
                        (d * sum + source + production[i] - uptake[i]) / (real * d + decay);
                    let after = ((1.0 - omega) * before + omega * balanced).max(0.0);
                    *out = after;
                    moved = moved.max((after - before).abs());
                }
                moved
            })
            .reduce(|| 0.0, f64::max);

        std::mem::swap(&mut self.values, &mut self.next);
        moved
    }

    /// What the held pixels hand the free ones per unit time, at the field as
    /// it stands: the diffusive flux across every face between a held pixel
    /// and a free one. At a settled field it is what the free pixels take up
    /// in all, less what they secrete.
    ///
    /// Summed in a fixed order, so that it is the same to the last bit on any
    /// number of threads.
    #[must_use]
    pub fn inflow(&self) -> f64 {
        let Some(span) = self.free_span() else {
            return 0.0;
        };
        let (w, h) = (self.width as usize, self.height as usize);
        let slice = w * h;
        let mut flux = 0.0;
        for z in span.z.0..=span.z.1 {
            for y in span.y.0..=span.y.1 {
                for x in span.x.0..=span.x.1 {
                    let i = z as usize * slice + y as usize * w + x as usize;
                    if self.held[i] {
                        continue;
                    }
                    let mut across = |j: usize| {
                        if self.held[j] {
                            flux += self.d * (self.values[j] - self.values[i]);
                        }
                    };
                    if x > 0 {
                        across(i - 1);
                    }
                    if x + 1 < self.width {
                        across(i + 1);
                    }
                    if y > 0 {
                        across(i - w);
                    }
                    if y + 1 < self.height {
                        across(i + w);
                    }
                    if z > 0 {
                        across(i - slice);
                    }
                    if z + 1 < self.depth {
                        across(i + slice);
                    }
                }
            }
        }
        flux
    }

    /// How many pixels are not held.
    #[must_use]
    pub fn free_count(&self) -> usize {
        self.held.iter().filter(|held| !**held).count()
    }

    /// Let go of every held pixel, before the sources are declared again.
    ///
    /// Holds accumulate, so that a bath and a vessel can coexist; the caller
    /// releases once per step and declares what holds this time.
    pub fn release_all(&mut self) {
        self.held.fill(false);
    }

    /// Hold every pixel `is_source` names at `value`, regardless of where it is.
    ///
    /// This is a source with a fixed place: a vessel, a needle, a bead of agar.
    /// Unlike [`ScalarField::hold_bath`] nothing has to reach it; it is held
    /// because it is there.
    pub fn hold_pixels(&mut self, value: f64, is_source: impl Fn(u32, u32) -> bool) {
        for y in 0..self.height {
            for x in 0..self.width {
                if is_source(x, y) {
                    let i = self.idx(x, y);
                    self.held[i] = true;
                    self.values[i] = value;
                }
            }
        }
    }

    /// Hold every pixel `is_source` names, by flat index, at `value`: as
    /// [`ScalarField::hold_pixels`], through every slice of a volume.
    pub fn hold_where(&mut self, value: f64, is_source: impl Fn(usize) -> bool) {
        for i in 0..self.values.len() {
            if is_source(i) {
                self.held[i] = true;
                self.values[i] = value;
            }
        }
    }

    /// Hold every pixel `level` gives a value at that value, each its own.
    ///
    /// A source whose strength varies from place to place: a vessel network
    /// that delivers more where more blood runs.
    pub fn hold_pixels_by(&mut self, level: impl Fn(u32, u32) -> Option<f64>) {
        for y in 0..self.height {
            for x in 0..self.width {
                if let Some(value) = level(x, y) {
                    let i = self.idx(x, y);
                    self.held[i] = true;
                    self.values[i] = value;
                }
            }
        }
    }

    /// Hold every pixel a bath can reach at `value`, on a flat field.
    ///
    /// `is_open` answers "can the bath occupy this pixel", which in practice means
    /// "is this pixel empty medium rather than tissue". Starting from the edge of
    /// the box, the bath spreads through open pixels and pins every one it reaches.
    ///
    /// Holds accumulate: call [`ScalarField::release_all`] first when the layout
    /// has changed, or pixels that were bath and are now tissue stay held.
    ///
    /// Reachability is what makes this right rather than merely convenient. A
    /// pocket of medium walled in by cells, which is what a resorbed necrotic core
    /// leaves behind, is never reached, so it does not turn into a source of
    /// oxygen in the middle of the tissue. Holding every empty pixel instead would
    /// re-supply precisely the region the model is meant to starve.
    ///
    /// The bath spreads through faces and not through corners, so a wall of cells
    /// touching only at their corners still seals a pocket. That is the
    /// conservative reading of an ambiguous case: leaving a pocket unsupplied
    /// costs a little realism, leaking a source into the tissue costs the result.
    pub fn hold_bath(&mut self, value: f64, is_open: impl Fn(u32, u32) -> bool) {
        let w = self.width as usize;
        self.hold_bath_at(value, |i| is_open((i % w) as u32, (i / w) as u32));
    }

    /// Hold every pixel a bath can reach at `value`, `is_open` asked by flat
    /// index: as [`ScalarField::hold_bath`], and in a volume from all six of
    /// its faces, spreading through the faces of voxels.
    ///
    /// Only the box around every closed pixel, one pixel wider, is flooded.
    /// Everything outside it is open, and reaches a face of the box along a
    /// straight line that never enters it, so it is held without being
    /// visited; the flood starts from the box's outer layer, open by
    /// construction and touching what is held. The pixels held are those a
    /// flood of the whole field would hold.
    pub fn hold_bath_at(&mut self, value: f64, is_open: impl Fn(usize) -> bool) {
        self.frontier.clear();
        let (w, h, d) = (
            self.width as usize,
            self.height as usize,
            self.depth as usize,
        );
        let slice = w * h;
        let at = |i: usize| [i % w, (i / w) % h, i / slice];

        let mut closed: Option<[(usize, usize); 3]> = None;
        for i in (0..self.values.len()).filter(|&i| !is_open(i)) {
            let p = at(i);
            let s = closed.get_or_insert([(p[0], p[0]), (p[1], p[1]), (p[2], p[2])]);
            for (axis, c) in s.iter_mut().zip(p) {
                *axis = (axis.0.min(c), axis.1.max(c));
            }
        }
        let limits = [w, h, d];
        let span = closed.map(|mut s| {
            for (axis, limit) in s.iter_mut().zip(limits) {
                *axis = (axis.0.saturating_sub(1), (axis.1 + 1).min(limit - 1));
            }
            s
        });
        let inside = |p: [usize; 3]| {
            span.is_some_and(|s| s.iter().zip(p).all(|(axis, c)| axis.0 <= c && c <= axis.1))
        };
        for i in 0..self.values.len() {
            if !inside(at(i)) {
                self.held[i] = true;
            }
        }
        if let Some(s) = span {
            // The outer layer of the box. A flat field's one slice is its
            // whole area, so its first and last slice are not a layer.
            let outer = |p: [usize; 3]| {
                p[0] == s[0].0
                    || p[0] == s[0].1
                    || p[1] == s[1].0
                    || p[1] == s[1].1
                    || (d > 1 && (p[2] == s[2].0 || p[2] == s[2].1))
            };
            for z in s[2].0..=s[2].1 {
                for y in s[1].0..=s[1].1 {
                    for x in s[0].0..=s[0].1 {
                        if outer([x, y, z]) {
                            self.reach(z * slice + y * w + x, &is_open);
                        }
                    }
                }
            }
        }

        while let Some(i) = self.frontier.pop() {
            let (x, y, z) = (i % w, (i / w) % h, i / slice);
            if x > 0 {
                self.reach(i - 1, &is_open);
            }
            if x + 1 < w {
                self.reach(i + 1, &is_open);
            }
            if y > 0 {
                self.reach(i - w, &is_open);
            }
            if y + 1 < h {
                self.reach(i + w, &is_open);
            }
            if z > 0 {
                self.reach(i - slice, &is_open);
            }
            if z + 1 < d {
                self.reach(i + slice, &is_open);
            }
        }

        for i in 0..self.values.len() {
            if self.held[i] {
                self.values[i] = value;
            }
        }
    }

    /// Mark one pixel as held by the bath, unless it is closed or already marked.
    fn reach(&mut self, i: usize, is_open: &impl Fn(usize) -> bool) {
        if self.held[i] || !is_open(i) {
            return;
        }
        self.held[i] = true;
        self.frontier.push(i);
    }

    fn idx(&self, x: u32, y: u32) -> usize {
        y as usize * self.width as usize + x as usize
    }

    /// Laplacian via a five point stencil, seven in a volume, with no-flux
    /// walls, at flat index `i`.
    ///
    /// At the edge of the box the missing neighbour is replaced by the pixel
    /// itself, which makes the flux through the wall exactly zero.
    fn laplacian(&self, i: usize) -> f64 {
        let (w, h) = (self.width as usize, self.height as usize);
        let slice = w * h;
        let (x, y, z) = (i % w, (i / w) % h, i / slice);
        let c = self.values[i];
        let at = |j: usize, there: bool| if there { self.values[j] } else { c };
        let mut lap = at(i.wrapping_sub(1), x > 0)
            + at(i + 1, x + 1 < w)
            + at(i.wrapping_sub(w), y > 0)
            + at(i + w, y + 1 < h)
            - 4.0 * c;
        if self.depth > 1 {
            lap += at(i.wrapping_sub(slice), z > 0) + at(i + slice, z + 1 < self.depth as usize)
                - 2.0 * c;
        }
        lap
    }

    /// Refuse a time step the explicit scheme cannot take.
    ///
    /// Past the limit the scheme does not fail, it oscillates and grows, which is
    /// the worst way for a simulation to be wrong: it still produces numbers.
    fn check_step(&self, dt: f64) -> Result<(), CellweaveError> {
        // Spelled out rather than negating a comparison, because NaN compares
        // false against everything and would otherwise slip past every check below.
        if dt.is_nan() || dt <= 0.0 {
            return Err(CellweaveError::Config(format!(
                "time step must be positive, got {dt}"
            )));
        }
        let limit = if self.depth > 1 {
            DIFFUSION_LIMIT * 2.0 / 3.0
        } else {
            DIFFUSION_LIMIT
        };
        if self.d * dt > limit {
            return Err(CellweaveError::Config(format!(
                "unstable time step: D·dt = {:.4} is above {limit} (D = {}, dt = {dt})",
                self.d * dt,
                self.d
            )));
        }
        if self.decay * dt > 1.0 {
            return Err(CellweaveError::Config(format!(
                "unstable time step: decay·dt = {:.4} is above 1 (decay = {}, dt = {dt})",
                self.decay * dt,
                self.decay
            )));
        }
        Ok(())
    }
}

impl DiffusionField for ScalarField {
    fn step(&mut self, dt: TimeStep) -> Result<(), CellweaveError> {
        let dt = dt.0;
        self.check_step(dt)?;

        for i in 0..self.values.len() {
            let c = self.values[i];
            if self.held[i] {
                // A bath pixel is a boundary condition, not an unknown.
                self.next[i] = c;
                continue;
            }
            let lap = self.laplacian(i);
            let change =
                self.d * lap - self.decay * c + self.source + self.production[i] - self.uptake[i];
            // The floor is physics, not a safety net: tissue cannot take up
            // what is not there. Diffusion and decay alone never go below
            // zero, the stability check sees to that.
            self.next[i] = (c + dt * change).max(0.0);
        }

        // Swap buffers, no allocation.
        std::mem::swap(&mut self.values, &mut self.next);
        Ok(())
    }

    fn concentration_at(&self, x: u32, y: u32) -> f64 {
        if x >= self.width || y >= self.height {
            return 0.0;
        }
        self.values[self.idx(x, y)]
    }

    fn dims(&self) -> (u32, u32) {
        (self.width, self.height)
    }
}

#[cfg(test)]
mod tests {
    #![allow(clippy::unwrap_used)]

    use super::*;

    /// Everything outside a disc of `radius` centred in a box of `size`.
    fn outside_disc(size: u32, radius: u32) -> impl Fn(u32, u32) -> bool {
        let centre = f64::from(size / 2);
        let r2 = f64::from(radius) * f64::from(radius);
        move |x, y| {
            let (dx, dy) = (f64::from(x) - centre, f64::from(y) - centre);
            dx * dx + dy * dy > r2
        }
    }

    /// Diffusion coefficient every disc test runs at.
    const D: f64 = 0.2;

    /// A consuming disc bathed on its surface, relaxed to its steady state.
    fn bathed_disc(size: u32, radius: u32, uptake: f64) -> ScalarField {
        let outside = outside_disc(size, radius);
        let mut field = ScalarField::new(size, size, D, 0.0, 1.0).unwrap();
        for y in 0..size {
            for x in 0..size {
                if !outside(x, y) {
                    field.set_uptake(x, y, uptake);
                }
            }
        }
        field.hold_bath(1.0, &outside);
        field.relax(TimeStep(1.0), 1e-11, 50_000).unwrap();
        field
    }

    /// Steady-state concentration at the middle of that disc.
    fn disc_centre(size: u32, radius: u32, uptake: f64) -> f64 {
        bathed_disc(size, radius, uptake).concentration_at(size / 2, size / 2)
    }

    #[test]
    fn the_size_of_the_box_does_not_change_the_tissue() {
        // The whole point of the bath. Supplying only the edge of the box turns
        // the medium into a funnel whose depth depends on how far the wall happens
        // to be, so doubling the grid changes the biology: at this uptake the same
        // disc reads 0.48 in a box of 30 and 0.22 in a box of 60. Held against the
        // surface of the tissue instead, the box stops being an actor.
        let small = disc_centre(30, 6, 0.004);
        let large = disc_centre(60, 6, 0.004);

        // Two grids mean two relaxation factors and two paths to the same
        // fixed point, so they agree to the solver's tolerance, not bit for bit.
        assert!(
            (small - large).abs() < 1e-9,
            "the box changed the answer: {small:.12} against {large:.12}, by {:e}",
            (small - large).abs()
        );
        // Both guards keep the comparison from being vacuous: a disc that never
        // depletes, or one flattened against zero, would agree for the wrong
        // reason.
        assert!(small < 0.9, "expected a depleted centre, got {small:.3}");
        assert!(
            small > 0.1,
            "expected the floor not to be reached, got {small:.3}"
        );
    }

    #[test]
    fn a_bathed_disc_matches_the_textbook_profile() {
        // A disc of radius R taking up q per unit time, bathed at c_s on its
        // surface, settles at c(r) = c_s - (q / 4D)(R² - r²), so the middle sits
        // a depth q·R²/(4D) below the surface. Reproducing a closed form nobody
        // tuned against is the strongest statement available about the solver.
        const RADIUS: u32 = 8;
        const UPTAKE: f64 = 0.003;

        let measured = disc_centre(32, RADIUS, UPTAKE);

        let depth = UPTAKE * f64::from(RADIUS * RADIUS) / (4.0 * D);
        let expected = 1.0 - depth;
        let error = (measured - expected).abs() / depth;
        // The disc is a staircase at this resolution, which costs a few percent
        // on where its surface really is; the error shrinks as the radius grows.
        assert!(
            error < 0.12,
            "measured {measured:.4} against {expected:.4}, {:.1}% off the drop",
            error * 100.0
        );
    }

    #[test]
    fn a_pocket_walled_in_by_tissue_is_never_reached() {
        let mut field = ScalarField::new(11, 11, 0.2, 0.0, 0.0).unwrap();
        // A square ring of tissue one pixel thick, five across.
        let wall = |x: u32, y: u32| {
            let (dx, dy) = (i64::from(x) - 5, i64::from(y) - 5);
            dx.abs().max(dy.abs()) == 3
        };

        field.hold_bath(1.0, |x, y| !wall(x, y));

        assert_eq!(field.concentration_at(5, 5), 0.0, "the pocket was supplied");
        assert_eq!(field.concentration_at(0, 0), 1.0, "the outside was not");
    }

    #[test]
    fn one_gap_in_the_wall_lets_the_bath_in() {
        let mut field = ScalarField::new(11, 11, 0.2, 0.0, 0.0).unwrap();
        let wall = |x: u32, y: u32| {
            let (dx, dy) = (i64::from(x) - 5, i64::from(y) - 5);
            dx.abs().max(dy.abs()) == 3 && !(x == 5 && y == 2)
        };

        field.hold_bath(1.0, |x, y| !wall(x, y));

        assert_eq!(field.concentration_at(5, 5), 1.0);
    }

    #[test]
    fn a_wall_that_only_touches_at_corners_still_seals() {
        // A diamond ring: its sides are diagonal chains of pixels meeting at their
        // corners only. Spreading through corners as well as faces would let the
        // bath slip between two of them and re-supply the pocket, so this test is
        // what pins the choice of four-connectivity.
        let mut field = ScalarField::new(11, 11, 0.2, 0.0, 0.0).unwrap();
        let wall = |x: u32, y: u32| {
            let (dx, dy) = (i64::from(x) - 5, i64::from(y) - 5);
            dx.abs() + dy.abs() == 3
        };

        field.hold_bath(1.0, |x, y| !wall(x, y));

        assert_eq!(
            field.concentration_at(5, 5),
            0.0,
            "the bath leaked diagonally"
        );
        assert_eq!(field.concentration_at(0, 0), 1.0);
    }

    #[test]
    fn a_closed_box_keeps_what_it_holds() {
        // No decay, no source, no bath, and walls that let nothing through: the
        // total must not move. It is the cheapest check that the stencil does not
        // quietly create or destroy matter at the edges.
        let mut field = ScalarField::new(20, 20, 0.2, 0.0, 0.0).unwrap();
        for y in 2..6u32 {
            for x in 2..6u32 {
                let i = field.idx(x, y);
                field.values[i] = 1.0;
            }
        }
        let before: f64 = field.values.iter().sum();

        for _ in 0..500 {
            field.step(TimeStep(1.0)).unwrap();
        }

        let after: f64 = field.values.iter().sum();
        assert!(
            (after - before).abs() < 1e-9,
            "total moved from {before:.9} to {after:.9}"
        );
    }

    #[test]
    fn an_unstable_time_step_is_refused() {
        let mut field = ScalarField::new(20, 20, 0.3, 0.0, 1.0).unwrap();

        let refused = field.step(TimeStep(1.0));

        assert!(matches!(refused, Err(CellweaveError::Config(_))));
    }

    #[test]
    fn a_decay_that_removes_more_than_everything_is_refused() {
        let mut field = ScalarField::new(20, 20, 0.1, 2.0, 1.0).unwrap();

        let refused = field.step(TimeStep(1.0));

        assert!(matches!(refused, Err(CellweaveError::Config(_))));
    }

    #[test]
    fn a_time_step_going_nowhere_is_refused() {
        let mut field = ScalarField::new(20, 20, 0.1, 0.0, 1.0).unwrap();

        assert!(field.step(TimeStep(0.0)).is_err());
        assert!(field.step(TimeStep(-1.0)).is_err());
    }

    #[test]
    fn a_field_with_no_pixels_is_refused() {
        assert!(ScalarField::new(0, 10, 0.2, 0.0, 1.0).is_err());
        assert!(ScalarField::new(10, 0, 0.2, 0.0, 1.0).is_err());
    }

    #[test]
    fn a_field_that_is_not_a_number_is_refused() {
        // NaN compares false against every bound, so a plain range check lets it
        // through and the whole grid turns into NaN on the first step.
        assert!(ScalarField::new(10, 10, f64::NAN, 0.0, 1.0).is_err());
        assert!(ScalarField::new(10, 10, 0.2, f64::NAN, 1.0).is_err());
        assert!(ScalarField::new(10, 10, f64::INFINITY, 0.0, 1.0).is_err());

        let mut field = ScalarField::new(10, 10, 0.2, 0.0, 1.0).unwrap();
        assert!(field.step(TimeStep(f64::NAN)).is_err());
    }

    #[test]
    fn field_decays_toward_zero_without_source() {
        let mut field = ScalarField::new(20, 20, 0.25, 0.1, 1.0).unwrap();

        for _ in 0..200 {
            field.step(TimeStep(0.1)).unwrap();
        }

        let centre = field.concentration_at(10, 10);
        assert!(centre < 0.5, "expected decay, centre = {centre:.3}");
    }

    #[test]
    fn concentration_never_negative() {
        // The stencil no longer clamps, so this is a check that refusing an
        // unstable step is enough on its own.
        let mut field = ScalarField::new(30, 30, 0.5, 2.0, 1.0).unwrap();

        for _ in 0..500 {
            field.step(TimeStep(0.05)).unwrap();
        }

        for y in 0..30u32 {
            for x in 0..30u32 {
                let c = field.concentration_at(x, y);
                assert!(c >= 0.0, "negative concentration at ({x},{y}): {c}");
            }
        }
    }

    #[test]
    fn the_number_of_threads_does_not_change_a_single_bit() {
        // A half sweep reads one grid and writes the other, so which thread
        // takes which rows cannot change a value. An odd-sized box, so the rows
        // do not split evenly between tasks.
        let solve = |threads| {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build()
                .unwrap()
                .install(|| bathed_disc(61, 20, 0.0005))
                .values
                .iter()
                .map(|v| v.to_bits())
                .collect::<Vec<u64>>()
        };
        let alone = solve(1);
        assert_eq!(alone, solve(3));
        assert_eq!(alone, solve(16));
    }

    #[test]
    fn what_the_bath_hands_a_settled_tissue_is_what_it_takes_up() {
        let field = bathed_disc(30, 6, 0.004);
        let taken = 0.004 * field.free_count() as f64;
        assert!(field.free_count() > 100);
        assert!(
            (field.inflow() - taken).abs() < 1e-6 * taken,
            "handed {} against {taken} taken up",
            field.inflow()
        );
    }

    #[test]
    fn relaxing_settles_the_field_and_says_how_long_it_took() {
        let mut field = bathed_disc(30, 6, 0.004);
        let settled: Vec<f64> = (0..900).map(|i| field.values[i]).collect();

        let more = field.relax(TimeStep(1.0), 1e-11, 50_000).unwrap();

        // Already settled: one step confirms it, and nothing moved.
        assert_eq!(more, 1);
        let drift = (0..900)
            .map(|i| (field.values[i] - settled[i]).abs())
            .fold(0.0_f64, f64::max);
        assert!(drift < 1e-6, "a settled field drifted by {drift:e}");
    }

    #[test]
    fn relaxing_refuses_to_pretend_when_it_cannot_settle_in_time() {
        let outside = outside_disc(30, 6);
        let mut field = ScalarField::new(30, 30, D, 0.0, 1.0).unwrap();
        for y in 0..30 {
            for x in 0..30 {
                if !outside(x, y) {
                    field.set_uptake(x, y, 0.004);
                }
            }
        }
        field.hold_bath(1.0, &outside);

        let refused = field.relax(TimeStep(1.0), 1e-9, 3);

        assert!(matches!(refused, Err(CellweaveError::Simulation(_))));
    }

    #[test]
    fn a_held_pixel_is_a_boundary_condition_and_does_not_move() {
        // The bath sits at 1.0 next to a pocket of tissue draining hard. Without
        // the hold, diffusion into the drain would pull the bath pixels down.
        let mut field = bathed_disc(30, 6, 0.05);

        for _ in 0..200 {
            field.step(TimeStep(1.0)).unwrap();
        }

        assert_eq!(field.concentration_at(0, 0), 1.0);
        assert_eq!(
            field.concentration_at(15, 3),
            1.0,
            "a bath pixel at the rim moved"
        );
        assert!(
            field.concentration_at(15, 15) < 0.5,
            "the drain should be starved"
        );
    }

    #[test]
    fn clearing_the_uptake_lets_the_field_refill() {
        let mut field = bathed_disc(30, 6, 0.01);
        assert!(field.concentration_at(15, 15) < 0.9);

        field.clear_uptake();
        field.relax(TimeStep(1.0), 1e-7, 50_000).unwrap();

        assert!(
            (field.concentration_at(15, 15) - 1.0).abs() < 1e-5,
            "with no sink and a bath at 1.0 the centre must return to 1.0"
        );
    }

    #[test]
    fn a_held_vessel_supplies_from_the_inside_out() {
        // No bath at the edge, a source in the middle, tissue draining all
        // around it: concentration must fall with distance from the source,
        // and the far corner, sealed from any bath, must be nearly empty.
        let mut field = ScalarField::new(60, 60, D, 0.0, 0.0).unwrap();
        let vessel = |x: u32, y: u32| {
            let (dx, dy) = (i64::from(x) - 30, i64::from(y) - 30);
            dx * dx + dy * dy <= 9
        };
        for y in 0..60 {
            for x in 0..60 {
                if !vessel(x, y) {
                    field.set_uptake(x, y, 0.002);
                }
            }
        }
        field.hold_pixels(1.0, vessel);
        field.relax(TimeStep(1.0), 1e-7, 50_000).unwrap();

        let at = |r: u32| field.concentration_at(30 + r, 30);
        assert_eq!(at(0), 1.0);
        assert!(
            at(5) > at(10) && at(10) > at(15),
            "not falling with distance"
        );
        assert!(
            field.concentration_at(1, 1) < 0.01,
            "the far corner should be starved"
        );
    }

    #[test]
    fn releasing_lets_a_former_bath_pixel_move_again() {
        let mut field = ScalarField::new(20, 20, D, 0.0, 1.0).unwrap();
        field.hold_bath(1.0, |_, _| true);
        for y in 0..20 {
            for x in 0..20 {
                field.set_uptake(x, y, 0.01);
            }
        }
        field.step(TimeStep(1.0)).unwrap();
        assert_eq!(
            field.concentration_at(10, 10),
            1.0,
            "held, so it must not move"
        );

        field.release_all();
        field.step(TimeStep(1.0)).unwrap();

        assert!(
            field.concentration_at(10, 10) < 1.0,
            "released, so the sink acts"
        );
    }

    #[test]
    fn a_bath_and_a_vessel_can_hold_at_once() {
        let mut field = ScalarField::new(30, 30, D, 0.0, 0.0).unwrap();
        let tissue = |x: u32, y: u32| (5..25).contains(&x) && (5..25).contains(&y);
        let vessel = |x: u32, y: u32| x == 15 && y == 15;

        field.hold_bath(1.0, |x, y| !tissue(x, y));
        field.hold_pixels(0.5, vessel);

        assert_eq!(field.concentration_at(0, 0), 1.0);
        assert_eq!(field.concentration_at(15, 15), 0.5);
        assert_eq!(
            field.concentration_at(10, 10),
            0.0,
            "tissue is held by neither"
        );
    }

    #[test]
    fn a_sweep_is_held_to_the_box_around_what_is_free() {
        let mut field = ScalarField::new(20, 20, D, 0.0, 0.0).unwrap();
        let tissue = |x: u32, y: u32| (4..9).contains(&x) && (11..14).contains(&y);
        field.hold_bath(1.0, |x, y| !tissue(x, y));
        assert_eq!(
            field.free_span(),
            Some(Span {
                x: (4, 8),
                y: (11, 13),
                z: (0, 0)
            })
        );

        field.hold_pixels(1.0, |_, _| true);
        assert_eq!(field.free_span(), None, "nothing left to balance");
        assert_eq!(field.relax(TimeStep(1.0), 1e-9, 10).unwrap(), 1);
    }

    #[test]
    fn the_box_swept_leaves_what_lies_outside_it_as_it_was_held() {
        // A starving disc in a bath: the sweep never touches the bath, which
        // must come out of the solve holding what it held going in.
        let mut field = ScalarField::new(40, 40, D, 0.0, 1.0).unwrap();
        let disc = |x: u32, y: u32| (x as i64 - 20).pow(2) + (y as i64 - 20).pow(2) <= 64;
        for y in 0..40 {
            for x in 0..40 {
                if disc(x, y) {
                    field.set_uptake(x, y, 0.001);
                }
            }
        }
        field.hold_bath(1.0, |x, y| !disc(x, y));
        field.relax(TimeStep(1.0), 1e-9, 100_000).unwrap();
        for y in 0..40 {
            for x in 0..40 {
                if !disc(x, y) {
                    assert_eq!(field.concentration_at(x, y), 1.0, "bath at ({x}, {y})");
                }
            }
        }
        assert!(field.concentration_at(20, 20) < 1.0, "the centre starves");
    }

    #[test]
    fn a_starving_ball_settles_where_the_sphere_formula_puts_it() {
        // Uptake u inside a ball of radius R in a bath held at one: the
        // centre sits at 1 - u R^2 / (6 D), where a disc's would be at
        // 1 - u R^2 / (4 D), a volume being fed from more directions.
        let (side, r, u) = (31_u32, 10.0_f64, 0.006);
        let mut field = ScalarField::with_depth(side, side, side, D, 0.0, 1.0).unwrap();
        let c = f64::from(side / 2);
        let inside = |i: usize| {
            let s = side as usize;
            let (x, y, z) = ((i % s) as f64, ((i / s) % s) as f64, (i / (s * s)) as f64);
            (x - c).powi(2) + (y - c).powi(2) + (z - c).powi(2) <= r * r
        };
        for i in 0..field.len() {
            if inside(i) {
                field.set_uptake_at(i, u);
            }
        }
        field.hold_bath_at(1.0, |i| !inside(i));
        field.relax(TimeStep(1.0), 1e-10, 1_000_000).unwrap();
        let centre = field
            .value_at((side as usize / 2) * (side as usize * side as usize + side as usize + 1));
        let expected = 1.0 - u * r * r / (6.0 * D);
        assert!(
            (centre - expected).abs() < 0.05,
            "centre {centre:.3}, the sphere formula {expected:.3}"
        );
    }

    #[test]
    fn a_bath_reaches_a_volume_through_its_six_faces() {
        // Open only along a column through the middle, which touches no side
        // of the box and both its front and back faces.
        let column = |i: usize| i % 64 == 4 * 8 + 4;
        let mut volume = ScalarField::with_depth(8, 8, 8, D, 0.0, 0.0).unwrap();
        volume.hold_bath_at(1.0, column);
        assert_eq!(
            volume.value_at(3 * 64 + 4 * 8 + 4),
            1.0,
            "reached from the faces"
        );
        assert_eq!(volume.value_at(3 * 64 + 4 * 8 + 3), 0.0, "closed");
        // In the plane the same pixel is walled in on every side.
        let mut flat = ScalarField::new(8, 8, D, 0.0, 0.0).unwrap();
        flat.hold_bath_at(1.0, column);
        assert_eq!(flat.value_at(4 * 8 + 4), 0.0, "a pocket in the plane");
    }
}
