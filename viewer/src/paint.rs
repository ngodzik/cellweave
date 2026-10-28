//! Turns the lattice into a picture: one texel block per pixel, each cell
//! filled with what it is doing and outlined where it meets something else.
//!
//! This is the colour of the surface under the lamp. What it gives off as
//! light is drawn separately, by the glow.

use bevy::asset::RenderAssetUsages;
use bevy::prelude::*;
use bevy::render::render_resource::{Extent3d, TextureDimension, TextureFormat};
use cellweave_core::{CellId, CloneId, DiffusionField, Seeding, SignalingNetwork};
use cellweave_engine::cpm::OBSTACLE;
use cellweave_engine::simulation::Running;
use cellweave_engine::{CpmLattice, Fate, Simulation};

use crate::runner::Runner;

/// Texels per lattice pixel along each side.
///
/// Enough that a membrane can be one texel wide, a quarter of a pixel, and
/// read as a line around a cell rather than as a cell of its own.
pub const SCALE: u32 = 4;

/// How bright the brightest emission is, in HDR units. Above one, it blooms.
pub const GLOW: f32 = 4.0;

/// What the window shows behind the tissue.
pub const BACKGROUND: Color = Color::srgb(0.004, 0.006, 0.014);

/// A cell alive and dividing.
pub const CYCLING: LinearRgba = LinearRgba::rgb(0.02, 0.62, 0.48);
/// A cell alive but arrested by hypoxia.
pub const QUIESCENT: LinearRgba = LinearRgba::rgb(0.90, 0.32, 0.02);
/// A dead cell, still giving up its pixels.
pub const NECROTIC: LinearRgba = LinearRgba::rgb(0.08, 0.05, 0.10);
/// A cell committed to differentiate, on its way out of the layer: pale, a
/// cell that has stopped being a progenitor.
pub const DIFFERENTIATING: LinearRgba = LinearRgba::rgb(0.55, 0.45, 0.28);
/// A cell in apoptosis, breaking up: violet, the colour of nothing else on
/// the canvas, so that where the stress is killing reads at a glance.
pub const APOPTOTIC: LinearRgba = LinearRgba::rgb(0.50, 0.14, 0.62);
/// The vessel wall, held at full oxygen. The brightest thing in the picture,
/// because everything else lives off it.
/// A lesion's cells in the clone view: a pale gold no clone's hue reaches,
/// so an early tumour stands out from the patchwork of clones around it.
pub const LESION: LinearRgba = LinearRgba::rgb(1.0, 0.82, 0.45);

pub const VESSEL: LinearRgba = LinearRgba::rgb(1.0, 0.08, 0.04);
/// The medium with no oxygen in it, and with all of it.
const MEDIUM_EMPTY: LinearRgba = LinearRgba::rgb(0.001, 0.002, 0.008);
const MEDIUM_FULL: LinearRgba = LinearRgba::rgb(0.002, 0.012, 0.032);

/// How much darker a membrane is than the cell it bounds.
const MEMBRANE: f32 = 0.12;

/// What the cells are coloured by.
#[derive(Resource, Debug, Default, Clone, Copy, PartialEq, Eq)]
pub enum View {
    /// Cycling, quiescent or necrotic: the decision each cell's network made.
    #[default]
    Fate,
    /// The level of one protein, by its index in the network.
    Protein(usize),
    /// The oxygen each pixel sits in, cells included.
    Oxygen,
    /// The clone each cell belongs to: which mutation it descends from.
    Clone,
    /// The lethal stress each pixel sits in, cells included.
    Stress,
}

/// The texture one run's tissue is painted into, and the material that
/// shows it.
#[derive(Component)]
pub struct Canvas {
    /// Written by [`repaint`].
    pub image: Handle<Image>,
    /// Shows the image as the surface's colour.
    pub material: Handle<StandardMaterial>,
    /// The generation last painted, and the view it was painted in.
    pub seen: Option<(u64, View)>,
}

/// A blank canvas for a lattice of this size, kept on the CPU so it can be
/// repainted every step.
pub fn canvas(width: u32, height: u32) -> Image {
    Image::new_fill(
        Extent3d {
            width: width * SCALE,
            height: height * SCALE,
            depth_or_array_layers: 1,
        },
        TextureDimension::D2,
        &[0, 0, 0, 255],
        TextureFormat::Rgba8UnormSrgb,
        RenderAssetUsages::MAIN_WORLD | RenderAssetUsages::RENDER_WORLD,
    )
}

/// Repaint each run's canvas when its simulation moved or the view changed.
///
/// A view changed while a simulation is away in the task pool is painted
/// when it comes back.
pub fn repaint(
    mut sides: Query<(&Runner, &mut Canvas)>,
    view: Res<View>,
    mut images: ResMut<Assets<Image>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    for (runner, mut canvas) in &mut sides {
        let Some(sim) = runner.idle() else { continue };
        let now = (runner.generation(), *view);
        if canvas.seen == Some(now) {
            continue;
        }
        let Some(mut image) = images.get_mut(&canvas.image) else {
            continue;
        };
        let Some(data) = image.data.as_mut() else {
            continue;
        };
        paint(sim, *view, data);
        // An image changed under a material does not always reach the
        // material's bind group; touching the material makes it pick the new
        // texture up.
        let _ = materials.get_mut(&canvas.material);
        canvas.seen = Some(now);
    }
}

/// Paint the whole lattice into `out`, an sRGB RGBA buffer [`SCALE`] times
/// the lattice on each side.
fn paint(sim: &Simulation<Running>, view: View, out: &mut [u8]) {
    let lattice = sim.lattice();
    let oxygen = sim.oxygen();
    let stress = |x: u32, y: u32| sim.stress().map_or(0.0, |f| f.concentration_at(x, y)) as f32;
    let (w, h) = lattice.dims();

    // Each cell's colour, worked out once rather than once per pixel.
    let cells: Vec<[u8; 4]> = (0..=lattice.cell_count() as u32)
        .map(|id| bytes(cell_colour(sim, CellId(id), view)))
        .collect();
    let vessel = bytes(VESSEL);

    for y in 0..h {
        for x in 0..w {
            let id = lattice.occupant(x, y);
            let fill = match id {
                OBSTACLE => vessel,
                _ if view == View::Stress => bytes(ramp(stress(x, y))),
                0 => bytes(MEDIUM_EMPTY.mix(&MEDIUM_FULL, oxygen.concentration_at(x, y) as f32)),
                _ if view == View::Oxygen => bytes(ramp(oxygen.concentration_at(x, y) as f32)),
                _ => cells.get(id as usize).copied().unwrap_or(vessel),
            };
            block(out, w, (x, y), fill, darken(fill), membrane(lattice, x, y));
        }
    }
}

/// Which sides of pixel `(x, y)` are a membrane: left, right, top, bottom.
///
/// Only a cell has membranes. A neighbour off the edge of the box counts as
/// the same cell, so the box itself is not outlined.
pub fn membrane(lattice: &CpmLattice, x: u32, y: u32) -> [bool; 4] {
    let (w, h) = lattice.dims();
    let id = lattice.occupant(x, y);
    if id == 0 || id == OBSTACLE {
        return [false; 4];
    }
    let differs = |nx: i64, ny: i64| {
        nx >= 0
            && ny >= 0
            && (nx as u32) < w
            && (ny as u32) < h
            && lattice.occupant(nx as u32, ny as u32) != id
    };
    let (xi, yi) = (i64::from(x), i64::from(y));
    [
        differs(xi - 1, yi),
        differs(xi + 1, yi),
        differs(xi, yi - 1),
        differs(xi, yi + 1),
    ]
}

/// Fill the [`SCALE`] by [`SCALE`] block of texels for pixel `at` of a
/// lattice `width` pixels wide, with `rim` on the sides `edges` marks.
pub fn block(
    out: &mut [u8],
    width: u32,
    at: (u32, u32),
    fill: [u8; 4],
    rim: [u8; 4],
    edges: [bool; 4],
) {
    let row = (width * SCALE * 4) as usize;
    let (x, y) = at;
    for sy in 0..SCALE {
        let start = (y * SCALE + sy) as usize * row + (x * SCALE * 4) as usize;
        for sx in 0..SCALE {
            let on_edge = (sx == 0 && edges[0])
                || (sx == SCALE - 1 && edges[1])
                || (sy == 0 && edges[2])
                || (sy == SCALE - 1 && edges[3]);
            let at = start + (sx * 4) as usize;
            out[at..at + 4].copy_from_slice(if on_edge { &rim } else { &fill });
        }
    }
}

/// The colour a cell is filled with, before its membrane is drawn.
pub fn cell_colour(sim: &Simulation<Running>, id: CellId, view: View) -> LinearRgba {
    let Some(fate) = sim.fate(id) else {
        return NECROTIC;
    };
    // A small, fixed variation per cell, so that neighbours in the same state
    // still read as separate cells and not as one coloured blob.
    let jitter = 0.88 + 0.24 * ((f64::from(id.0) * 0.618_034).fract() as f32);
    let base = match (view, fate) {
        (_, Fate::Necrotic) => return NECROTIC,
        (View::Protein(i), _) => {
            let level = sim
                .network(id)
                .and_then(|net| net.state().values.get(i))
                .map_or(0.0, |c| c.0 as f32);
            ramp(level)
        }
        (View::Clone, fate) => {
            let seeding = sim.scenario().seeding;
            let colour = if sim.in_lesion(id) {
                LESION
            } else {
                sim.clone_of(id)
                    .map_or(NECROTIC, |c| clone_colour(c, seeding))
            };
            // A clone keeps its hue whatever its cells are doing; an arrested
            // or departing cell is a dimmer shade of it, so a clone that
            // escapes arrest is also a clone that lights up.
            if matches!(
                fate,
                Fate::Quiescent | Fate::Differentiating | Fate::Apoptotic
            ) {
                colour * 0.35
            } else {
                colour
            }
        }
        (_, Fate::Cycling) => CYCLING,
        (_, Fate::Quiescent) => QUIESCENT,
        (_, Fate::Differentiating) => DIFFERENTIATING,
        (_, Fate::Apoptotic) => APOPTOTIC,
    };
    base * jitter
}

/// The colour of a clone. Grown from one cell, the founding clone is a cold,
/// muted blue, so that every mutant stands out against it. In a sheet every
/// cell founds a clone and none is special. Other clones step round the hue
/// circle by the golden angle, which keeps clones founded one after the other
/// far apart.
pub fn clone_colour(clone: CloneId, seeding: Seeding) -> LinearRgba {
    if clone.0 == 0 && seeding == Seeding::Single {
        return LinearRgba::rgb(0.10, 0.20, 0.42);
    }
    let hue = (f64::from(clone.0) * 137.507_764).rem_euclid(360.0) as f32;
    Color::hsl(hue, 0.85, 0.45).to_linear()
}

/// A dark-to-bright scale for a value between zero and one: violet, red,
/// amber, then near white, so that a high level is also a bright one.
pub fn ramp(value: f32) -> LinearRgba {
    const STOPS: [(f32, LinearRgba); 5] = [
        (0.0, LinearRgba::rgb(0.01, 0.005, 0.04)),
        (0.3, LinearRgba::rgb(0.20, 0.02, 0.35)),
        (0.6, LinearRgba::rgb(0.85, 0.10, 0.16)),
        (0.85, LinearRgba::rgb(1.0, 0.52, 0.05)),
        (1.0, LinearRgba::rgb(1.0, 0.95, 0.60)),
    ];
    let v = value.clamp(0.0, 1.0);
    STOPS
        .windows(2)
        .find(|pair| v <= pair[1].0)
        .map_or(STOPS[4].1, |pair| {
            let (a, b) = (pair[0], pair[1]);
            a.1.mix(&b.1, (v - a.0) / (b.0 - a.0))
        })
}

/// A linear colour as the sRGB bytes the canvas holds.
pub fn bytes(colour: LinearRgba) -> [u8; 4] {
    Srgba::from(colour).to_u8_array()
}

fn darken(fill: [u8; 4]) -> [u8; 4] {
    let d = |c: u8| (f32::from(c) * MEMBRANE) as u8;
    [d(fill[0]), d(fill[1]), d(fill[2]), 255]
}
