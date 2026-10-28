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
            // A cell leaving the layer keeps its grip on its neighbours until
            // it is gone, so it costs what it cost while it was dividing.
            CellKind::Normal | CellKind::Differentiated => 2,
            // A cell breaking up holds on to its neighbours no better than a
            // dead one.
            CellKind::Necrotic | CellKind::Apoptotic => 3,
        }
    }
}

/// Below this many pixels a cell is too small to be worth splitting.
const MIN_DIVIDE_VOLUME: u32 = 8;

/// Identifier of pixels no cell can ever take: a vessel wall, a piece of
/// glass. Not a cell, so it has no record, no energy, no volume; a flip into or
/// out of it is simply never proposed.
pub const OBSTACLE: u32 = u32::MAX;

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
                let free = self.grid[[row as usize, col as usize]] != OBSTACLE;
                if dx * dx + dy * dy <= r2 && free {
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

    /// Add a rectangular cell of `kind`, `width` by `height` pixels, with its
    /// top left corner at `corner`, clipped to the box and to free pixels.
    ///
    /// Rectangles tile without gaps or overlap, which is how a sheet of tissue
    /// fills the box at the start. Returns the new `CellId`.
    pub fn add_block(&mut self, kind: CellKind, corner: Pos2, width: u32, height: u32) -> CellId {
        let id = self.cells.len() as u32;
        let mut vol = 0u32;
        for row in corner.y..(corner.y + height).min(self.height) {
            for col in corner.x..(corner.x + width).min(self.width) {
                let pixel = &mut self.grid[[row as usize, col as usize]];
                if *pixel == 0 {
                    *pixel = id;
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

    /// Fill a disc with obstacle, overwriting whatever was there.
    ///
    /// Any cell pixels overwritten are taken off their cells' volumes, so the
    /// records stay true to the grid. Called before cells are placed for a
    /// vessel the tissue grows around, and afterwards for one that grows
    /// through the tissue.
    pub fn add_obstacle_disc(&mut self, centre: Pos2, radius: u32) {
        let cx = i64::from(centre.x);
        let cy = i64::from(centre.y);
        let r2 = i64::from(radius * radius);
        // Only the square around the disc can be in it.
        let rows = centre.y.saturating_sub(radius)..(centre.y + radius + 1).min(self.height);
        let cols = centre.x.saturating_sub(radius)..(centre.x + radius + 1).min(self.width);
        for row in rows {
            for col in cols.clone() {
                let dx = i64::from(col) - cx;
                let dy = i64::from(row) - cy;
                if dx * dx + dy * dy > r2 {
                    continue;
                }
                let was = self.grid[[row as usize, col as usize]];
                if was != 0 && was != OBSTACLE {
                    self.cells[was as usize].volume.0 -= 1;
                }
                self.grid[[row as usize, col as usize]] = OBSTACLE;
            }
        }
    }

    /// Fill a disc with obstacle, pushing the cells in it aside rather than
    /// taking their pixels.
    ///
    /// The disc's cell pixels are shoved one at a time, the nearest to the
    /// centre first. Each moves one pixel outward, along the line from the
    /// centre through it, and so does every cell pixel beyond it on that
    /// line, up to the first pixel of medium, which the last of them fills.
    /// No cell changes size. The centre itself is shoved along `heading`, the
    /// way the disc is moving. A pixel whose line meets an obstacle or the
    /// edge of the box before any medium cannot be shoved, and is taken off
    /// its cell as [`CpmLattice::add_obstacle_disc`] takes it. Returns how
    /// many pixels were taken.
    pub fn push_obstacle_disc(&mut self, centre: Pos2, radius: u32, heading: (f64, f64)) -> u32 {
        let (cx, cy) = (i64::from(centre.x), i64::from(centre.y));
        let r2 = i64::from(radius * radius);
        let rows = centre.y.saturating_sub(radius)..(centre.y + radius + 1).min(self.height);
        let cols = centre.x.saturating_sub(radius)..(centre.x + radius + 1).min(self.width);
        let mut disc: Vec<(i64, i64)> = rows
            .flat_map(|y| cols.clone().map(move |x| (i64::from(x), i64::from(y))))
            .filter(|&(x, y)| (x - cx).pow(2) + (y - cy).pow(2) <= r2)
            .collect();
        // Stable, so that pixels as near as each other go in reading order.
        disc.sort_by_key(|&(x, y)| (x - cx).pow(2) + (y - cy).pow(2));
        let mut taken = 0;
        for (x, y) in disc {
            let id = self.grid[[y as usize, x as usize]];
            if id == OBSTACLE {
                continue;
            }
            let direction = if (x, y) == (cx, cy) {
                heading
            } else {
                ((x - cx) as f64, (y - cy) as f64)
            };
            if id != 0 && !self.shove((x, y), direction) {
                self.cells[id as usize].volume.0 -= 1;
                taken += 1;
            }
            self.grid[[y as usize, x as usize]] = OBSTACLE;
        }
        taken
    }

    /// Move the cell pixel at `from`, and every cell pixel beyond it along
    /// `direction`, one pixel further, into the first pixel of medium on the
    /// way. The line steps one pixel at a time along its longer axis. Leaves
    /// `from` holding what it held, about to be overwritten, and returns
    /// whether medium was found before an obstacle or the edge of the box.
    fn shove(&mut self, from: (i64, i64), direction: (f64, f64)) -> bool {
        let longer = direction.0.abs().max(direction.1.abs());
        if !(longer.is_finite() && longer > 0.0) {
            return false;
        }
        let (sx, sy) = (direction.0 / longer, direction.1 / longer);
        let at = |k: i64| {
            (
                from.0 + (k as f64 * sx).round() as i64,
                from.1 + (k as f64 * sy).round() as i64,
            )
        };
        let inside = |(x, y): (i64, i64)| {
            x >= 0 && y >= 0 && x < i64::from(self.width) && y < i64::from(self.height)
        };
        let mut end = 1;
        loop {
            let pixel = at(end);
            if !inside(pixel) {
                return false;
            }
            match self.grid[[pixel.1 as usize, pixel.0 as usize]] {
                0 => break,
                OBSTACLE => return false,
                _ => end += 1,
            }
        }
        for k in (1..=end).rev() {
            let (to, behind) = (at(k), at(k - 1));
            self.grid[[to.1 as usize, to.0 as usize]] =
                self.grid[[behind.1 as usize, behind.0 as usize]];
        }
        true
    }

    /// Give an obstacle pixel back to the medium, where a vessel was taken
    /// back, for cells to move into. Any other pixel is left as it is.
    pub fn clear_obstacle(&mut self, x: u32, y: u32) {
        if let Some(pixel) = self.grid.get_mut([y as usize, x as usize])
            && *pixel == OBSTACLE
        {
            *pixel = 0;
        }
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

        if src == dst || src == OBSTACLE || dst == OBSTACLE {
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
        if id == OBSTACLE {
            // No record. A boundary against an obstacle costs what a boundary
            // against medium costs, so an obstacle is a wall a cell can rest on.
            return CellKind::Medium;
        }
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

    /// Biological kind of a cell.
    pub fn cell_kind(&self, id: CellId) -> Option<CellKind> {
        self.cells.get(id.0 as usize).map(|r| r.kind)
    }

    /// Change a cell's kind, for instance when it dies.
    ///
    /// The medium and unknown identifiers are ignored.
    pub fn set_kind(&mut self, id: CellId, kind: CellKind) {
        if id.0 == 0 {
            return;
        }
        if let Some(r) = self.cells.get_mut(id.0 as usize) {
            r.kind = kind;
        }
    }

    /// The identifier occupying one pixel: 0 for medium, [`OBSTACLE`] for a
    /// pixel no cell can take.
    ///
    /// Coordinates outside the grid read as medium.
    pub fn occupant(&self, x: u32, y: u32) -> u32 {
        if x >= self.width || y >= self.height {
            return 0;
        }
        self.grid[[y as usize, x as usize]]
    }

    /// Width and height.
    pub fn dims(&self) -> (u32, u32) {
        (self.width, self.height)
    }

    /// Apply the volume parameters a cell's signaling network asks for.
    ///
    /// This is the upward half of the coupling: protein state decides how large a
    /// cell tries to be. Adhesion is still a global table, so the contact fields
    /// of [`cellweave_core::CellMechanics`] have no effect yet.
    ///
    /// An unknown identifier is ignored rather than reported, since the caller
    /// iterates over cells it just read from this lattice.
    pub fn set_volume_target(&mut self, id: CellId, target: Volume, lambda: f64) {
        if let Some(r) = self.cells.get_mut(id.0 as usize) {
            r.target_volume = target;
            r.lambda_volume = lambda;
        }
    }

    /// Split a cell in two along a straight line through its centre of mass.
    ///
    /// Pixels on one side of the line keep the parent identity, the rest go to a
    /// new cell of the same kind that inherits the parent's volume parameters.
    /// `angle` orients the cut, in radians.
    ///
    /// Returns the daughter's identifier, or `None` when the cell is the medium,
    /// unknown, too small to split, or when the line leaves one side empty, in
    /// which case nothing is changed.
    ///
    /// This scans the whole lattice once, so it is meant to run every so often
    /// rather than inside the Monte Carlo loop.
    pub fn divide(&mut self, id: CellId, angle: f64) -> Option<CellId> {
        if id.0 == 0 {
            return None;
        }
        let parent = self.cells.get(id.0 as usize)?.clone();

        let mut pixels = Vec::with_capacity(parent.volume.0 as usize);
        for y in 0..self.height {
            for x in 0..self.width {
                if self.grid[[y as usize, x as usize]] == id.0 {
                    pixels.push((x, y));
                }
            }
        }
        // The grid is the truth about volume; the record is a cache of it.
        let total = pixels.len() as u32;
        if total < MIN_DIVIDE_VOLUME {
            return None;
        }

        let n = f64::from(total);
        let cx = pixels.iter().map(|(x, _)| f64::from(*x)).sum::<f64>() / n;
        let cy = pixels.iter().map(|(_, y)| f64::from(*y)).sum::<f64>() / n;
        let (nx, ny) = (angle.cos(), angle.sin());

        let daughter_id = self.cells.len() as u32;
        let mut moved = 0u32;
        for (x, y) in &pixels {
            if (f64::from(*x) - cx) * nx + (f64::from(*y) - cy) * ny > 0.0 {
                self.grid[[*y as usize, *x as usize]] = daughter_id;
                moved += 1;
            }
        }

        // A line that leaves one side empty is not a division. Put it back.
        if moved == 0 || moved == total {
            for (x, y) in &pixels {
                self.grid[[*y as usize, *x as usize]] = id.0;
            }
            return None;
        }

        self.cells[id.0 as usize].volume = Volume(total - moved);
        self.cells.push(CellRecord {
            kind: parent.kind,
            volume: Volume(moved),
            target_volume: parent.target_volume,
            lambda_volume: parent.lambda_volume,
        });
        Some(CellId(daughter_id))
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

    /// Pixel count per identifier, read from the grid rather than the records.
    /// Obstacles have no record and are skipped.
    fn recount(lat: &CpmLattice) -> Vec<u32> {
        let mut counts = vec![0u32; lat.cells.len()];
        for id in &lat.grid {
            if *id != OBSTACLE {
                counts[*id as usize] += 1;
            }
        }
        counts
    }

    /// Every record's volume matches what the grid actually holds.
    fn records_match_grid(lat: &CpmLattice) -> bool {
        recount(lat)
            .iter()
            .zip(&lat.cells)
            .skip(1)
            .all(|(counted, record)| *counted == record.volume.0)
    }

    #[test]
    fn a_cell_grows_toward_a_raised_target() {
        let mut lat = CpmLattice::new(120, 120, 12.0).unwrap();
        let id = lat.add_cell(CellKind::Tumor, Pos2 { x: 60, y: 60 }, 12);
        let start = lat.volume(id).unwrap().0;

        lat.set_volume_target(id, Volume(start * 2), 20.0);
        let mut rng = StdRng::seed_from_u64(7);
        for _ in 0..400 {
            lat.monte_carlo_step(&mut rng);
        }

        let grown = lat.volume(id).unwrap().0;
        assert!(
            grown > start + start / 4,
            "expected clear growth, went from {start} to {grown}"
        );
        assert!(records_match_grid(&lat));
    }

    #[test]
    fn dividing_hands_every_pixel_to_exactly_one_of_the_two() {
        let mut lat = CpmLattice::new(60, 60, 10.0).unwrap();
        let parent = lat.add_cell(CellKind::Tumor, Pos2 { x: 30, y: 30 }, 10);
        let before = lat.volume(parent).unwrap().0;

        let daughter = lat.divide(parent, 0.3).expect("a round cell splits");

        let (p, d) = (
            lat.volume(parent).unwrap().0,
            lat.volume(daughter).unwrap().0,
        );
        assert_eq!(p + d, before, "pixels were created or lost");
        assert!(p > 0 && d > 0);
        // A line through the centre of a disc splits it near the middle.
        assert!((p as i64 - d as i64).unsigned_abs() < before as u64 / 5);
        assert!(records_match_grid(&lat));
    }

    #[test]
    fn the_daughter_takes_the_parents_kind_and_targets() {
        let mut lat = CpmLattice::new(60, 60, 10.0).unwrap();
        let parent = lat.add_cell(CellKind::Normal, Pos2 { x: 30, y: 30 }, 10);
        lat.set_volume_target(parent, Volume(500), 33.0);

        let daughter = lat.divide(parent, 1.0).unwrap();

        assert_eq!(daughter, CellId(2));
        assert_eq!(lat.cell_kind(daughter), Some(CellKind::Normal));
        assert_eq!(lat.cells[2].target_volume, Volume(500));
        assert_eq!(lat.cells[2].lambda_volume, 33.0);
        assert_eq!(lat.cell_count(), 2);
    }

    #[test]
    fn a_cell_too_small_to_split_is_left_alone() {
        let mut lat = CpmLattice::new(30, 30, 10.0).unwrap();
        let tiny = lat.add_cell(CellKind::Tumor, Pos2 { x: 15, y: 15 }, 1);
        let grid_before = lat.grid.clone();

        assert_eq!(lat.divide(tiny, 0.5), None);

        assert_eq!(lat.grid, grid_before);
        assert_eq!(lat.cell_count(), 1);
    }

    #[test]
    fn the_medium_cannot_be_divided() {
        let mut lat = CpmLattice::new(30, 30, 10.0).unwrap();
        lat.add_cell(CellKind::Tumor, Pos2 { x: 15, y: 15 }, 4);
        let grid_before = lat.grid.clone();

        assert_eq!(lat.divide(CellId(0), 0.5), None);
        assert_eq!(lat.divide(CellId(42), 0.5), None);

        assert_eq!(lat.grid, grid_before);
    }

    #[test]
    fn a_cut_that_leaves_one_side_empty_is_rolled_back() {
        // A one pixel wide vertical bar, cut by a line whose normal points along x:
        // every pixel projects to exactly zero, so nothing lands on the far side.
        let mut lat = CpmLattice::new(20, 20, 10.0).unwrap();
        lat.cells.push(CellRecord {
            kind: CellKind::Tumor,
            volume: Volume(10),
            target_volume: Volume(10),
            lambda_volume: 1.0,
        });
        for y in 5..15u32 {
            lat.grid[[y as usize, 10]] = 1;
        }
        let grid_before = lat.grid.clone();

        assert_eq!(lat.divide(CellId(1), 0.0), None);

        assert_eq!(lat.grid, grid_before, "a refused cut must change nothing");
        assert_eq!(lat.cell_count(), 1);
        assert!(records_match_grid(&lat));
    }

    #[test]
    fn bookkeeping_survives_division_followed_by_stepping() {
        let mut lat = CpmLattice::new(80, 80, 10.0).unwrap();
        let a = lat.add_cell(CellKind::Tumor, Pos2 { x: 40, y: 40 }, 12);
        let mut rng = StdRng::seed_from_u64(3);

        let b = lat.divide(a, 0.7).unwrap();
        for _ in 0..200 {
            lat.monte_carlo_step(&mut rng);
        }
        lat.divide(b, 2.1);
        for _ in 0..200 {
            lat.monte_carlo_step(&mut rng);
        }

        assert!(records_match_grid(&lat), "records drifted from the grid");
    }

    #[test]
    fn an_obstacle_is_never_taken_and_never_moves() {
        let mut lat = CpmLattice::new(60, 60, 30.0).unwrap();
        lat.add_obstacle_disc(Pos2 { x: 30, y: 30 }, 5);
        let cell = lat.add_cell(CellKind::Tumor, Pos2 { x: 30, y: 30 }, 12);
        let obstacle_before: Vec<(usize, usize)> = (0..60)
            .flat_map(|y| (0..60).map(move |x| (x, y)))
            .filter(|&(x, y)| lat.grid[[y, x]] == OBSTACLE)
            .collect();
        assert!(!obstacle_before.is_empty());

        // Hot and hungry: a high temperature and a target far above the volume
        // are what would push a cell into anything it could take.
        lat.set_volume_target(cell, Volume(2000), 50.0);
        let mut rng = StdRng::seed_from_u64(11);
        for _ in 0..300 {
            lat.monte_carlo_step(&mut rng);
        }

        for (x, y) in &obstacle_before {
            assert_eq!(lat.grid[[*y, *x]], OBSTACLE, "obstacle lost at ({x}, {y})");
        }
        assert!(records_match_grid(&lat));
    }

    #[test]
    fn a_cell_placed_over_an_obstacle_wraps_around_it() {
        let mut lat = CpmLattice::new(40, 40, 10.0).unwrap();
        lat.add_obstacle_disc(Pos2 { x: 20, y: 20 }, 4);
        let obstacle = lat.grid.iter().filter(|&&v| v == OBSTACLE).count() as u32;

        let cell = lat.add_cell(CellKind::Tumor, Pos2 { x: 20, y: 20 }, 8);

        let taken = lat.grid.iter().filter(|&&v| v != 0).count() as u32;
        assert_eq!(lat.volume(cell).unwrap().0, taken - obstacle);
        assert!(records_match_grid(&lat));
    }

    #[test]
    fn a_cell_aiming_at_nothing_shrinks_to_nothing() {
        // What resorption relies on: a target of zero drains a cell pixel by
        // pixel until none is left, and the record says so.
        let mut lat = CpmLattice::new(50, 50, 10.0).unwrap();
        let cell = lat.add_cell(CellKind::Necrotic, Pos2 { x: 25, y: 25 }, 6);
        lat.set_volume_target(cell, Volume(0), 50.0);
        let mut rng = StdRng::seed_from_u64(5);

        for _ in 0..200 {
            lat.monte_carlo_step(&mut rng);
        }

        assert_eq!(lat.volume(cell), Some(Volume(0)));
        assert!(
            lat.grid.iter().all(|&v| v != cell.0),
            "pixels outlived the record"
        );
        assert!(records_match_grid(&lat));
    }

    #[test]
    fn a_gone_cell_stays_gone_while_another_grows_through_its_old_place() {
        let mut lat = CpmLattice::new(50, 50, 10.0).unwrap();
        let gone = lat.add_cell(CellKind::Necrotic, Pos2 { x: 20, y: 25 }, 5);
        let alive = lat.add_cell(CellKind::Tumor, Pos2 { x: 32, y: 25 }, 5);
        lat.set_volume_target(gone, Volume(0), 50.0);
        lat.set_volume_target(alive, Volume(400), 50.0);
        let mut rng = StdRng::seed_from_u64(9);

        for _ in 0..400 {
            lat.monte_carlo_step(&mut rng);
        }

        assert_eq!(lat.volume(gone), Some(Volume(0)));
        assert!(
            lat.volume(alive).unwrap().0 > 200,
            "the living cell should have spread"
        );
        assert!(records_match_grid(&lat));
    }

    #[test]
    fn a_disc_pushed_into_tissue_moves_cells_aside_and_keeps_every_size() {
        // Three cells side by side, with medium all around them.
        let mut lat = CpmLattice::new(60, 60, 10.0).unwrap();
        let ids: Vec<CellId> = (0..3)
            .map(|k| {
                lat.add_block(
                    CellKind::Tumor,
                    Pos2 {
                        x: 15 + 10 * k,
                        y: 25,
                    },
                    10,
                    10,
                )
            })
            .collect();
        let sizes = |lat: &CpmLattice| -> Vec<u32> {
            ids.iter().map(|&id| lat.volume(id).unwrap().0).collect()
        };
        let before = sizes(&lat);
        let medium = |lat: &CpmLattice| lat.grid.iter().filter(|&&p| p == 0).count();
        let free = medium(&lat);

        let taken = lat.push_obstacle_disc(Pos2 { x: 30, y: 30 }, 2, (1.0, 0.0));
        assert_eq!(taken, 0, "there was room to push into");
        assert_eq!(sizes(&lat), before, "no cell lost a pixel");
        assert!(records_match_grid(&lat));
        let disc = lat.grid.iter().filter(|&&p| p == OBSTACLE).count();
        assert_eq!(disc, 13, "a disc of radius two");
        assert_eq!(medium(&lat), free - 13, "what was pushed filled the medium");
    }

    #[test]
    fn a_pixel_with_nothing_free_beyond_it_is_taken() {
        // One cell filling the box: nowhere to push it.
        let mut lat = CpmLattice::new(20, 20, 10.0).unwrap();
        let id = lat.add_block(CellKind::Tumor, Pos2 { x: 0, y: 0 }, 20, 20);
        let taken = lat.push_obstacle_disc(Pos2 { x: 10, y: 10 }, 2, (1.0, 0.0));
        assert_eq!(taken, 13);
        assert_eq!(lat.volume(id).unwrap().0, 400 - 13);
        assert!(records_match_grid(&lat));
    }

    #[test]
    fn a_push_stops_at_an_obstacle_and_takes_the_pixel_instead() {
        // One cell from x = 3 to 19, a wall at x = 20, medium on both sides
        // of the pair. Pushed from a disc at (16, 10): every line to the
        // right meets the wall, the vertical and diagonal ones the edge of
        // the box, and only the two leftward along the row reach medium.
        let mut lat = CpmLattice::new(40, 20, 10.0).unwrap();
        for y in 0..20 {
            lat.add_obstacle_disc(Pos2 { x: 20, y }, 0);
        }
        let id = lat.add_block(CellKind::Tumor, Pos2 { x: 3, y: 0 }, 17, 20);
        let left = |lat: &CpmLattice| {
            (0..20)
                .flat_map(|y| (0..3).map(move |x| (x, y)))
                .filter(|&(x, y)| lat.occupant(x, y) == 0)
                .count()
        };
        let free = left(&lat);
        let taken = lat.push_obstacle_disc(Pos2 { x: 16, y: 10 }, 2, (1.0, 0.0));
        assert_eq!(taken, 11);
        assert_eq!(lat.volume(id).unwrap().0, 17 * 20 - 11);
        assert_eq!(left(&lat), free - 2, "the two pushed went left");
        assert!(records_match_grid(&lat));
    }
}
