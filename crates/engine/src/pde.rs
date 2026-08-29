//! Scalar extracellular diffusion field, explicit finite differences.
//!
//! Solves ∂c/∂t = D·∇²c − λ·c + S on a 2D grid.

use cellweave_core::{CellweaveError, DiffusionField, TimeStep};

/// Scalar 2D diffusion field.
///
/// Boundary condition: Dirichlet (constant value at edges), representing
/// a well-perfused vascular boundary.
pub struct ScalarField {
    width: u32,
    height: u32,
    /// Current concentration on each pixel.
    values: Vec<f64>,
    /// Scratch buffer for the next time step.
    next: Vec<f64>,
    /// Diffusion coefficient (pixels²/s).
    d: f64,
    /// First-order decay rate (1/s).
    decay: f64,
    /// Dirichlet boundary value (kept for potential future write-back).
    #[allow(dead_code)]
    boundary: f64,
    /// Uniform source term (mol/L/s).
    source: f64,
}

impl ScalarField {
    /// Create a field initialised to `initial_value` everywhere.
    pub fn new(width: u32, height: u32, d: f64, decay: f64, boundary: f64) -> Self {
        let n = (width * height) as usize;
        Self {
            width,
            height,
            values: vec![boundary; n],
            next: vec![boundary; n],
            d,
            decay,
            boundary,
            source: 0.0,
        }
    }

    /// Set a spatially uniform source term.
    pub fn with_source(mut self, source: f64) -> Self {
        self.source = source;
        self
    }

    fn idx(&self, x: u32, y: u32) -> usize {
        y as usize * self.width as usize + x as usize
    }

    /// Laplacian at interior pixel (x, y) via 5-point stencil.
    fn laplacian(&self, x: u32, y: u32) -> f64 {
        let c = self.values[self.idx(x, y)];
        let left = self.values[self.idx(x - 1, y)];
        let right = self.values[self.idx(x + 1, y)];
        let up = self.values[self.idx(x, y - 1)];
        let down = self.values[self.idx(x, y + 1)];
        left + right + up + down - 4.0 * c
    }
}

impl DiffusionField for ScalarField {
    fn step(&mut self, dt: TimeStep) -> Result<(), CellweaveError> {
        let dt = dt.0;
        let d = self.d;

        // Update interior pixels; boundary stays at Dirichlet value.
        for y in 1..(self.height - 1) {
            for x in 1..(self.width - 1) {
                let i = self.idx(x, y);
                let c = self.values[i];
                let lap = self.laplacian(x, y);
                self.next[i] = c + dt * (d * lap - self.decay * c + self.source);
                // Prevent negative concentrations.
                self.next[i] = self.next[i].max(0.0);
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
    use super::*;

    #[test]
    fn field_decays_toward_zero_without_source() {
        let mut field = ScalarField::new(20, 20, 0.25, 0.1, 0.0);
        // Set interior to 1.0 manually.
        for y in 1..19u32 {
            for x in 1..19u32 {
                let i = field.idx(x, y);
                field.values[i] = 1.0;
            }
        }
        for _ in 0..200 {
            field.step(TimeStep(0.1)).unwrap();
        }
        // Interior should have decayed significantly.
        let centre = field.concentration_at(10, 10);
        assert!(centre < 0.5, "expected decay, centre = {:.3}", centre);
    }

    #[test]
    fn boundary_remains_constant() {
        let mut field = ScalarField::new(20, 20, 0.25, 0.0, 1.0);
        for _ in 0..100 {
            field.step(TimeStep(0.1)).unwrap();
        }
        assert_eq!(field.concentration_at(0, 0), 1.0);
        assert_eq!(field.concentration_at(19, 19), 1.0);
    }

    #[test]
    fn concentration_never_negative() {
        let mut field = ScalarField::new(30, 30, 0.5, 2.0, 0.0);
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
}
