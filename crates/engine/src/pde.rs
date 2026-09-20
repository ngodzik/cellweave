//! Scalar extracellular diffusion field, explicit finite differences.
//!
//! Solves ∂c/∂t = D·∇²c − λ·c + S on a 2D grid.
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

/// Scalar 2D diffusion field with no-flux walls.
pub struct ScalarField {
    width: u32,
    height: u32,
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
    /// Scratch: breadth-first frontier for the bath fill.
    frontier: Vec<(u32, u32)>,
}

/// Largest `D · dt` an explicit five point stencil stays stable at, on a unit grid.
const DIFFUSION_LIMIT: f64 = 0.25;

impl ScalarField {
    /// Create a field holding `initial` everywhere.
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
        if width == 0 || height == 0 {
            return Err(CellweaveError::Config(format!(
                "field dimensions must be positive, got {width}x{height}"
            )));
        }
        if !d.is_finite() || !decay.is_finite() || d < 0.0 || decay < 0.0 {
            return Err(CellweaveError::Config(format!(
                "diffusion and decay must be finite and not negative, got D = {d}, decay = {decay}"
            )));
        }
        let n = width as usize * height as usize;
        Ok(Self {
            width,
            height,
            values: vec![initial; n],
            next: vec![initial; n],
            d,
            decay,
            source: 0.0,
            held: vec![false; n],
            uptake: vec![0.0; n],
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

    /// Forget every sink, before the layout of the tissue is read in again.
    pub fn clear_uptake(&mut self) {
        self.uptake.fill(0.0);
    }

    /// Iterate until the field stops moving, and say how many steps it took.
    ///
    /// Converged means no pixel moved by more than `tolerance` over the last
    /// step. Started from the previous layout's solution this takes a few dozen
    /// steps; started cold it can take thousands, which is why `max_steps` is a
    /// hard limit that turns a field which will not settle into an error rather
    /// than a silent partial answer.
    ///
    /// # Errors
    ///
    /// Whatever [`DiffusionField::step`] refuses, and [`CellweaveError::Simulation`]
    /// when `max_steps` pass without convergence.
    pub fn relax(
        &mut self,
        dt: TimeStep,
        tolerance: f64,
        max_steps: usize,
    ) -> Result<usize, CellweaveError> {
        for done in 1..=max_steps {
            self.step(dt)?;
            // After the swap, `next` holds the previous values.
            let moved = self
                .values
                .iter()
                .zip(&self.next)
                .map(|(now, before)| (now - before).abs())
                .fold(0.0_f64, f64::max);
            if moved < tolerance {
                return Ok(done);
            }
        }
        Err(CellweaveError::Simulation(format!(
            "diffusion did not settle within {max_steps} steps at tolerance {tolerance}"
        )))
    }

    /// Hold every pixel a bath can reach at `value`.
    ///
    /// `is_open` answers "can the bath occupy this pixel", which in practice means
    /// "is this pixel empty medium rather than tissue". Starting from the edge of
    /// the box, the bath spreads through open pixels and pins every one it reaches.
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
        self.held.fill(false);
        self.frontier.clear();

        let (last_x, last_y) = (self.width - 1, self.height - 1);
        for x in 0..self.width {
            self.reach(x, 0, &is_open);
            self.reach(x, last_y, &is_open);
        }
        for y in 0..self.height {
            self.reach(0, y, &is_open);
            self.reach(last_x, y, &is_open);
        }

        while let Some((x, y)) = self.frontier.pop() {
            if x > 0 {
                self.reach(x - 1, y, &is_open);
            }
            if x < last_x {
                self.reach(x + 1, y, &is_open);
            }
            if y > 0 {
                self.reach(x, y - 1, &is_open);
            }
            if y < last_y {
                self.reach(x, y + 1, &is_open);
            }
        }

        for i in 0..self.values.len() {
            if self.held[i] {
                self.values[i] = value;
            }
        }
    }

    /// Mark one pixel as held by the bath, unless it is closed or already marked.
    fn reach(&mut self, x: u32, y: u32, is_open: &impl Fn(u32, u32) -> bool) {
        let i = self.idx(x, y);
        if self.held[i] || !is_open(x, y) {
            return;
        }
        self.held[i] = true;
        self.frontier.push((x, y));
    }

    fn idx(&self, x: u32, y: u32) -> usize {
        y as usize * self.width as usize + x as usize
    }

    /// Laplacian via a five point stencil, with no-flux walls.
    ///
    /// At the edge of the box the missing neighbour is replaced by the pixel
    /// itself, which makes the flux through the wall exactly zero.
    fn laplacian(&self, x: u32, y: u32) -> f64 {
        let c = self.values[self.idx(x, y)];
        let left = self.values[self.idx(x.saturating_sub(1), y)];
        let right = self.values[self.idx((x + 1).min(self.width - 1), y)];
        let up = self.values[self.idx(x, y.saturating_sub(1))];
        let down = self.values[self.idx(x, (y + 1).min(self.height - 1))];
        left + right + up + down - 4.0 * c
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
        if self.d * dt > DIFFUSION_LIMIT {
            return Err(CellweaveError::Config(format!(
                "unstable time step: D·dt = {:.4} is above {DIFFUSION_LIMIT} (D = {}, dt = {dt})",
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

        for y in 0..self.height {
            for x in 0..self.width {
                let i = self.idx(x, y);
                let c = self.values[i];
                if self.held[i] {
                    // A bath pixel is a boundary condition, not an unknown.
                    self.next[i] = c;
                    continue;
                }
                let lap = self.laplacian(x, y);
                let change = self.d * lap - self.decay * c + self.source - self.uptake[i];
                // The floor is physics, not a safety net: tissue cannot take up
                // what is not there. Diffusion and decay alone never go below
                // zero, the stability check sees to that.
                self.next[i] = (c + dt * change).max(0.0);
            }
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
        field.relax(TimeStep(1.0), 1e-7, 50_000).unwrap();
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

        assert!(
            (small - large).abs() < 1e-9,
            "the box changed the answer: {small:.6} against {large:.6}"
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
    fn relaxing_settles_the_field_and_says_how_long_it_took() {
        let mut field = bathed_disc(30, 6, 0.004);
        let settled: Vec<f64> = (0..900).map(|i| field.values[i]).collect();

        let more = field.relax(TimeStep(1.0), 1e-7, 50_000).unwrap();

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
}
