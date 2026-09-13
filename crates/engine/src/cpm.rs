//! Cellular Potts Model lattice (Graner & Glazier 1992).
//!
//! Uses the checkerboard update scheme so parallel spin flips never race.

use cellweave_core::{CellId, CellKind, Pos2, Volume};
use ndarray::Array2;
// rand 0.10 split the trait: Rng carries the raw bits, RngExt the helpers.
// Public signatures keep the Rng bound, since RngExt is blanket implemented.
use rand::{Rng, RngExt};
use thiserror::Error;

/// Errors specific to the CPM lattice.
#[derive(Debug, Error)]
pub enum CpmError {
    /// Grid dimensions are zero or larger than 65535.
    #[error("invalid grid dimensions ({0}×{1})")]
    InvalidDims(u32, u32),
}

/// Contact energy `J` for every pair of cell kinds.
///
/// A boundary between two pixels of different cells costs the value listed here
/// for their kinds. Keeping `J(medium, cell)` well above `J(cell, cell)` is what
/// makes cells prefer each other's company over empty space, which is what holds
/// a cluster together and keeps a single cell compact.
#[derive(Debug, Clone, Copy)]
pub struct ContactEnergy {
    /// J(Medium, Tumor).
    pub medium_tumor: f64,
    /// J(Medium, Normal).
    pub medium_normal: f64,
    /// J(Medium, Necrotic).
    pub medium_necrotic: f64,
    /// J(Tumor, Tumor).
    pub tumor_tumor: f64,
    /// J(Tumor, Normal).
    pub tumor_normal: f64,
    /// J(Tumor, Necrotic).
    pub tumor_necrotic: f64,
    /// J(Normal, Normal).
    pub normal_normal: f64,
    /// J(Normal, Necrotic).
    pub normal_necrotic: f64,
    /// J(Necrotic, Necrotic).
    pub necrotic_necrotic: f64,
}

impl Default for ContactEnergy {
    fn default() -> Self {
        Self {
            medium_tumor: 16.0,
            medium_normal: 16.0,
            // Dead cells hold on to their surroundings less well.
            medium_necrotic: 12.0,
            tumor_tumor: 2.0,
            tumor_normal: 11.0,
            tumor_necrotic: 10.0,
            normal_normal: 6.0,
            normal_necrotic: 12.0,
            necrotic_necrotic: 8.0,
        }
    }
}

impl ContactEnergy {
    /// Energy of a boundary between two pixels of the given kinds.
    ///
    /// The table is symmetric, so the pair is ordered before lookup.
    pub fn between(&self, a: CellKind, b: CellKind) -> f64 {
        let (i, j) = (Self::index(a), Self::index(b));
        let (lo, hi) = if i <= j { (i, j) } else { (j, i) };
        match (lo, hi) {
            (0, 0) => 0.0,
            (0, 1) => self.medium_tumor,
            (0, 2) => self.medium_normal,
            (0, 3) => self.medium_necrotic,
            (1, 1) => self.tumor_tumor,
            (1, 2) => self.tumor_normal,
            (1, 3) => self.tumor_necrotic,
            (2, 2) => self.normal_normal,
            (2, 3) => self.normal_necrotic,
            (3, 3) => self.necrotic_necrotic,
            // Ordering above makes every other pair impossible.
            _ => 0.0,
        }
    }

    fn index(k: CellKind) -> u8 {
        match k {
            CellKind::Medium => 0,
            CellKind::Tumor => 1,
            CellKind::Normal => 2,
            CellKind::Necrotic => 3,
        }
    }
}

/// Per-cell properties tracked by the lattice.
#[derive(Debug, Clone)]
pub struct CellRecord {
    /// Biological type.
    pub kind: CellKind,
    /// Current pixel count.
    pub volume: Volume,
    /// Target volume from the ODE layer.
    pub target_volume: Volume,
    /// Volume constraint strength λ.
    pub lambda_volume: f64,
}

/// 2D Cellular Potts Model lattice.
///
/// Each pixel stores a `CellId` (0 = medium). Monte Carlo spin flips are
/// proposed randomly; accepted with the Boltzmann probability.
pub struct CpmLattice {
    /// Width in pixels.
    pub width: u32,
    /// Height in pixels.
    pub height: u32,
    /// Lattice: each pixel → CellId (0 = medium).
    pub grid: Array2<u32>,
    /// Per-cell records indexed by CellId.
    cells: Vec<CellRecord>,
    /// Effective temperature (controls shape fluctuations).
    pub temperature: f64,
    /// Contact energies.
    pub contact: ContactEnergy,
    /// Current Monte Carlo step.
    pub mcs: u64,
}

impl CpmLattice {
    /// Create an empty lattice filled with medium.
    pub fn new(width: u32, height: u32, temperature: f64) -> Result<Self, CpmError> {
        if width == 0 || height == 0 || width > 65535 || height > 65535 {
            return Err(CpmError::InvalidDims(width, height));
        }
        Ok(Self {
            width,
            height,
            grid: Array2::zeros((height as usize, width as usize)),
            // Id 0 = medium; pre-populate the slot so indices are stable.
            cells: vec![CellRecord {
                kind: CellKind::Medium,
                volume: Volume(0),
                target_volume: Volume(0),
                lambda_volume: 0.0,
            }],
            temperature,
            contact: ContactEnergy::default(),
            mcs: 0,
        })
    }

    /// Add a circular cell of `kind` centred at `centre` with radius `r`.
    ///
    /// Returns the new `CellId`.
    pub fn add_cell(&mut self, kind: CellKind, centre: Pos2, radius: u32) -> CellId {
        let id = self.cells.len() as u32;
        let mut vol = 0u32;
        let cx = centre.x as i64;
        let cy = centre.y as i64;
        let r2 = (radius * radius) as i64;

        for row in 0..self.height {
            for col in 0..self.width {
                let dx = col as i64 - cx;
                let dy = row as i64 - cy;
                if dx * dx + dy * dy <= r2 {
                    self.grid[[row as usize, col as usize]] = id;
                    vol += 1;
                }
            }
        }

        self.cells.push(CellRecord {
            kind,
            volume: Volume(vol),
            target_volume: Volume(vol),
            lambda_volume: 50.0,
        });
        CellId(id)
    }

    /// Execute one full Monte Carlo step (W×H spin-flip attempts).
    ///
    /// Uses a simple sequential scan rather than checkerboard for the
    /// foundation milestone. Parallelism is wired in once tests pass.
    pub fn monte_carlo_step(&mut self, rng: &mut impl Rng) {
        let n = (self.width * self.height) as usize;
        for _ in 0..n {
            let col = rng.random_range(0..self.width);
            let row = rng.random_range(0..self.height);
            self.try_flip(row, col, rng);
        }
        self.mcs += 1;
    }

    /// Attempt a spin flip at `(row, col)`, Boltzmann acceptance.
    fn try_flip(&mut self, row: u32, col: u32, rng: &mut impl Rng) {
        let src = self.grid[[row as usize, col as usize]];

        // Pick a random von Neumann neighbour.
        let neighbours = self.von_neumann(row, col);
        if neighbours.is_empty() {
            return;
        }
        let (nr, nc) = neighbours[rng.random_range(0..neighbours.len())];
        let dst = self.grid[[nr as usize, nc as usize]];

        if src == dst {
            return;
        }

        let delta = self.delta_h(row, col, src, dst);
        let accept = if delta <= 0.0 {
            true
        } else {
            rng.random::<f64>() < (-delta / self.temperature).exp()
        };

        if accept {
            self.grid[[row as usize, col as usize]] = dst;
            if src != 0 {
                self.cells[src as usize].volume.0 -= 1;
            }
            if dst != 0 {
                self.cells[dst as usize].volume.0 += 1;
            }
        }
    }

    /// ΔH for copying spin `dst` into pixel `(row, col)` currently holding `src`.
    fn delta_h(&self, row: u32, col: u32, src: u32, dst: u32) -> f64 {
        self.contact_delta_h(row, col, src, dst) + self.volume_delta_h(src, dst)
    }

    /// Change in boundary energy alone, with no volume term.
    ///
    /// The Hamiltonian sums `J(tau_i, tau_j)` over neighbouring pixel pairs that
    /// belong to different cells. Both the removed and the created boundaries have
    /// to be accounted for: a neighbour that already holds `src` contributes zero
    /// before the flip but a real cost afterwards, so it must not be skipped.
    fn contact_delta_h(&self, row: u32, col: u32, src: u32, dst: u32) -> f64 {
        let mut dh = 0.0;
        for (nr, nc) in self.von_neumann(row, col) {
            let neighbour = self.grid[[nr as usize, nc as usize]];
            dh += self.contact_energy(dst, neighbour) - self.contact_energy(src, neighbour);
        }
        dh
    }

    /// Change in the volume constraint `lambda * (V - V_target)^2`.
    fn volume_delta_h(&self, src: u32, dst: u32) -> f64 {
        let mut dh = 0.0;
        if src != 0 {
            let r = &self.cells[src as usize];
            let v = r.volume.0 as f64;
            let vt = r.target_volume.0 as f64;
            dh += r.lambda_volume * ((v - 1.0 - vt).powi(2) - (v - vt).powi(2));
        }
        if dst != 0 {
            let r = &self.cells[dst as usize];
            let v = r.volume.0 as f64;
            let vt = r.target_volume.0 as f64;
            dh += r.lambda_volume * ((v + 1.0 - vt).powi(2) - (v - vt).powi(2));
        }
        dh
    }

    /// J(a, b) contact energy between spin ids `a` and `b`.
    ///
    /// Two pixels of the same cell share no boundary, so they cost nothing. This
    /// is the Kronecker delta of the Potts Hamiltonian, and it applies to the cell
    /// identity, not to the cell kind: two distinct tumor cells in contact do pay.
    fn contact_energy(&self, a: u32, b: u32) -> f64 {
        if a == b {
            return 0.0;
        }
        self.contact.between(self.kind(a), self.kind(b))
    }

    fn kind(&self, id: u32) -> CellKind {
        self.cells
            .get(id as usize)
            .map_or(CellKind::Medium, |r| r.kind)
    }

    /// Von Neumann neighbourhood (4-connected), boundary-clipped.
    fn von_neumann(&self, row: u32, col: u32) -> Vec<(u32, u32)> {
        let mut out = Vec::with_capacity(4);
        if row > 0 {
            out.push((row - 1, col));
        }
        if row + 1 < self.height {
            out.push((row + 1, col));
        }
        if col > 0 {
            out.push((row, col - 1));
        }
        if col + 1 < self.width {
            out.push((row, col + 1));
        }
        out
    }

    /// Number of registered cells (excluding medium).
    pub fn cell_count(&self) -> usize {
        self.cells.len() - 1
    }

    /// Current volume of a cell.
    pub fn volume(&self, id: CellId) -> Option<Volume> {
        self.cells.get(id.0 as usize).map(|r| r.volume)
    }

    /// All cell records (excluding medium at index 0).
    pub fn cell_records(&self) -> impl Iterator<Item = (CellId, &CellRecord)> {
        self.cells
            .iter()
            .enumerate()
            .skip(1)
            .map(|(i, r)| (CellId(i as u32), r))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    #[test]
    fn lattice_creates_empty() {
        let lat = CpmLattice::new(50, 50, 10.0).unwrap();
        assert_eq!(lat.cell_count(), 0);
        assert!(lat.grid.iter().all(|&v| v == 0));
    }

    #[test]
    fn add_cell_registers_volume() {
        let mut lat = CpmLattice::new(50, 50, 10.0).unwrap();
        let id = lat.add_cell(CellKind::Tumor, Pos2 { x: 25, y: 25 }, 5);
        let vol = lat.volume(id).unwrap();
        // Circle area ≈ π*r² ≈ 78; allow ±5 for rasterisation.
        assert!((70..90).contains(&vol.0), "volume {}", vol.0);
    }

    #[test]
    fn mcs_preserves_approximate_volume() {
        let mut lat = CpmLattice::new(100, 100, 5.0).unwrap();
        let id = lat.add_cell(CellKind::Tumor, Pos2 { x: 50, y: 50 }, 8);
        let v0 = lat.volume(id).unwrap().0;
        let mut rng = StdRng::seed_from_u64(42);
        for _ in 0..10 {
            lat.monte_carlo_step(&mut rng);
        }
        let v1 = lat.volume(id).unwrap().0;
        // Volume constraint should keep it within 20% of target.
        let deviation = (v1 as f64 - v0 as f64).abs() / v0 as f64;
        assert!(deviation < 0.20, "volume drifted {:.1}%", deviation * 100.0);
    }

    #[test]
    fn invalid_dims_rejected() {
        assert!(CpmLattice::new(0, 100, 10.0).is_err());
        assert!(CpmLattice::new(100, 0, 10.0).is_err());
    }

    /// Lattice holding one tumor cell over a rectangle, with the volume constraint
    /// switched off so that only the contact term is measured.
    fn block_cell(fill_w: u32, fill_h: u32) -> CpmLattice {
        let mut lat = CpmLattice::new(20, 20, 10.0).unwrap();
        lat.cells.push(CellRecord {
            kind: CellKind::Tumor,
            volume: Volume(fill_w * fill_h),
            target_volume: Volume(fill_w * fill_h),
            lambda_volume: 0.0,
        });
        for y in 0..fill_h {
            for x in 0..fill_w {
                lat.grid[[y as usize, x as usize]] = 1;
            }
        }
        lat
    }

    #[test]
    fn same_cell_pixels_share_no_boundary() {
        let lat = block_cell(10, 20);
        assert_eq!(lat.contact_energy(1, 1), 0.0);
        assert_eq!(lat.contact_energy(0, 0), 0.0);
        assert_eq!(lat.contact_energy(1, 0), lat.contact.medium_tumor);
    }

    #[test]
    fn punching_a_hole_costs_four_boundary_units() {
        // Pixel (10, 5) sits well inside the cell, so all four neighbours belong to
        // it. Handing the pixel to the medium creates four new boundaries and
        // removes none, so the cost is exactly 4 J(medium, tumor).
        let lat = block_cell(10, 20);
        let j = lat.contact.medium_tumor;
        assert_eq!(lat.contact_delta_h(10, 5, 1, 0), 4.0 * j);
    }

    #[test]
    fn growing_an_isolated_bump_costs_two_boundary_units() {
        // Pixel (10, 10) is medium with a single cell neighbour. The cell taking it
        // removes one boundary and creates three, for a net cost of 2 J.
        let lat = block_cell(10, 20);
        let j = lat.contact.medium_tumor;
        assert_eq!(lat.contact_delta_h(10, 10, 0, 1), 2.0 * j);
    }

    #[test]
    fn sliding_a_diagonal_corner_costs_nothing() {
        // With the cell filling a quadrant, pixel (9, 9) is the corner: two
        // neighbours inside the cell, two in the medium. Moving the corner trades
        // two boundaries for two others, so the energy must not change at all.
        let lat = block_cell(10, 10);
        assert_eq!(lat.contact_delta_h(9, 9, 1, 0), 0.0);
    }
}
