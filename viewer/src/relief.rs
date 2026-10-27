//! The tissue as a surface in relief: every cell a dome, every membrane a
//! groove, lit from the side so that the shape of each cell reads as shadow
//! and not only as colour.
//!
//! Heights are worked out when a batch comes back from the pool, and eased
//! towards every frame, so the tissue moves continuously even though the
//! simulation only advances a few times a second.

use bevy::asset::RenderAssetUsages;
use bevy::mesh::{Indices, PrimitiveTopology};
use bevy::prelude::*;
use cellweave_core::CellId;
use cellweave_engine::cpm::OBSTACLE;
use cellweave_engine::simulation::Running;
use cellweave_engine::{Fate, Simulation};

use crate::runner::Runner;

/// How far from its membrane a cell reaches its full height, in pixels.
/// About half a cell's radius, so a cell is a dome and not a mesa.
const SHOULDER: f32 = 5.0;

/// Full height of a cell in each fate, and of the vessel, in pixels.
///
/// A cycling cell stands tallest; an arrested one a little lower; a dead one
/// has collapsed. The vessel stands above them all, as the thing they live
/// off.
const CYCLING_HEIGHT: f32 = 5.0;
const QUIESCENT_HEIGHT: f32 = 3.8;
const NECROTIC_HEIGHT: f32 = 1.0;
/// A cell leaving the layer rises out of it, the way a differentiating cell
/// in an epithelium moves up out of the basal layer.
const DIFFERENTIATING_HEIGHT: f32 = 8.0;
const VESSEL_HEIGHT: f32 = 9.0;

/// How quickly the surface catches up with the simulation, per second. At
/// eight, a change is nine tenths done in under a third of a second.
const EASE: f32 = 8.0;

/// The surface, and where it is heading.
#[derive(Resource)]
pub struct Relief {
    /// The mesh being reshaped.
    pub mesh: Handle<Mesh>,
    width: usize,
    height: usize,
    /// Height of each pixel in the last state drawn from the simulation.
    target: Vec<f32>,
    /// Height of each pixel on screen.
    current: Vec<f32>,
    /// One when the relief is shown, zero when it is flattened; eased to.
    lift: f32,
    /// Whether the relief is shown.
    pub raised: bool,
    /// The generation `target` was drawn from.
    seen: Option<u64>,
}

impl Relief {
    /// A flat surface for a lattice this size, raised once the first state
    /// arrives.
    pub fn new(mesh: Handle<Mesh>, width: u32, height: u32) -> Self {
        let n = (width * height) as usize;
        Self {
            mesh,
            width: width as usize,
            height: height as usize,
            target: vec![0.0; n],
            current: vec![0.0; n],
            lift: 0.0,
            raised: true,
            seen: None,
        }
    }
}

/// A flat grid with one vertex at the centre of every pixel, the texture
/// spread across it so that each pixel's texels sit around its vertex.
///
/// Kept on the CPU as well as the GPU, since it is rewritten every frame.
pub fn surface(width: u32, height: u32) -> Mesh {
    let (w, h) = (width as usize, height as usize);
    let mut positions = Vec::with_capacity(w * h);
    let mut uvs = Vec::with_capacity(w * h);
    for y in 0..h {
        for x in 0..w {
            positions.push(at(x, y, w, h, 0.0));
            uvs.push([(x as f32 + 0.5) / w as f32, (y as f32 + 0.5) / h as f32]);
        }
    }
    let mut indices = Vec::with_capacity((w - 1) * (h - 1) * 6);
    for y in 0..h - 1 {
        for x in 0..w - 1 {
            let a = (y * w + x) as u32;
            let (b, c, d) = (a + 1, a + w as u32, a + w as u32 + 1);
            // Counter-clockwise seen from above, so the faces point up.
            indices.extend_from_slice(&[a, c, b, b, c, d]);
        }
    }
    Mesh::new(
        PrimitiveTopology::TriangleList,
        RenderAssetUsages::MAIN_WORLD | RenderAssetUsages::RENDER_WORLD,
    )
    .with_inserted_attribute(Mesh::ATTRIBUTE_POSITION, positions)
    .with_inserted_attribute(Mesh::ATTRIBUTE_NORMAL, vec![[0.0, 1.0, 0.0]; w * h])
    .with_inserted_attribute(Mesh::ATTRIBUTE_UV_0, uvs)
    .with_inserted_indices(Indices::U32(indices))
}

/// Where pixel `(x, y)` sits in the world, at height `lift`: centred on the
/// origin, rows running along +Z.
fn at(x: usize, y: usize, w: usize, h: usize, lift: f32) -> [f32; 3] {
    [
        x as f32 + 0.5 - w as f32 / 2.0,
        lift,
        y as f32 + 0.5 - h as f32 / 2.0,
    ]
}

/// Work out the heights of the latest state, when there is a new one.
pub fn retarget(runner: Res<Runner>, relief: Option<ResMut<Relief>>) {
    let Some(mut relief) = relief else { return };
    let Some(sim) = runner.idle() else { return };
    if relief.seen == Some(runner.generation()) {
        return;
    }
    relief.seen = Some(runner.generation());
    relief.target = heights(sim);
}

/// The height of every pixel: zero in the medium, and inside a cell or the
/// vessel a dome, rising from its boundary over [`SHOULDER`] pixels to the
/// full height of what it is.
fn heights(sim: &Simulation<Running>) -> Vec<f32> {
    let lattice = sim.lattice();
    let (w, h) = lattice.dims();
    let (w, h) = (w as usize, h as usize);
    let id = |x: usize, y: usize| lattice.occupant(x as u32, y as u32);

    // Distance from each pixel to the nearest pixel that is not the same
    // occupant, by two chamfer passes: down and right, then up and left.
    let mut distance = vec![f32::MAX; w * h];
    for y in 0..h {
        for x in 0..w {
            let here = id(x, y);
            if here == 0 {
                distance[y * w + x] = 0.0;
                continue;
            }
            let boundary = (x == 0 || id(x - 1, y) != here)
                || (x + 1 == w || id(x + 1, y) != here)
                || (y == 0 || id(x, y - 1) != here)
                || (y + 1 == h || id(x, y + 1) != here);
            if boundary {
                distance[y * w + x] = 0.0;
            }
        }
    }
    const DIAGONAL: f32 = std::f32::consts::SQRT_2;
    let forward = [
        (-1, 0, 1.0),
        (0, -1, 1.0),
        (-1, -1, DIAGONAL),
        (1, -1, DIAGONAL),
    ];
    let backward = [
        (1, 0, 1.0),
        (0, 1, 1.0),
        (1, 1, DIAGONAL),
        (-1, 1, DIAGONAL),
    ];
    let mut pass = |order: &mut dyn Iterator<Item = (usize, usize)>, steps: &[(i64, i64, f32)]| {
        for (x, y) in order {
            let mut best = distance[y * w + x];
            for &(dx, dy, cost) in steps {
                let (nx, ny) = (x as i64 + dx, y as i64 + dy);
                if nx >= 0 && ny >= 0 && (nx as usize) < w && (ny as usize) < h {
                    best = best.min(distance[ny as usize * w + nx as usize] + cost);
                }
            }
            distance[y * w + x] = best;
        }
    };
    pass(
        &mut (0..h).flat_map(|y| (0..w).map(move |x| (x, y))),
        &forward,
    );
    pass(
        &mut (0..h).rev().flat_map(|y| (0..w).rev().map(move |x| (x, y))),
        &backward,
    );

    let crest = |cell: u32| match cell {
        0 => 0.0,
        OBSTACLE => VESSEL_HEIGHT,
        _ => match sim.fate(CellId(cell)) {
            Some(Fate::Cycling) => CYCLING_HEIGHT,
            Some(Fate::Quiescent) => QUIESCENT_HEIGHT,
            Some(Fate::Differentiating) => DIFFERENTIATING_HEIGHT,
            Some(Fate::Necrotic) | None => NECROTIC_HEIGHT,
        },
    };
    let mut raw = vec![0.0; w * h];
    for y in 0..h {
        for x in 0..w {
            // A quarter circle in profile: steep at the membrane, flat on top.
            let t = (distance[y * w + x] / SHOULDER).min(1.0);
            let dome = (1.0 - (1.0 - t) * (1.0 - t)).sqrt();
            raw[y * w + x] = crest(id(x, y)) * dome;
        }
    }
    // A cell's outline is a staircase of pixels. Half of a three by three
    // average rounds the steps off; all of it would fill the grooves between
    // cells.
    let mut out = vec![0.0; w * h];
    for y in 0..h {
        for x in 0..w {
            let (x0, x1) = (x.saturating_sub(1), (x + 1).min(w - 1));
            let (y0, y1) = (y.saturating_sub(1), (y + 1).min(h - 1));
            let mut sum = 0.0;
            for yy in y0..=y1 {
                for xx in x0..=x1 {
                    sum += raw[yy * w + xx];
                }
            }
            let mean = sum / ((x1 - x0 + 1) * (y1 - y0 + 1)) as f32;
            out[y * w + x] = 0.5 * (raw[y * w + x] + mean);
        }
    }
    out
}

/// H raises or flattens the relief, unless egui is taking the keys.
pub fn toggle(
    keys: Res<ButtonInput<KeyCode>>,
    relief: Option<ResMut<Relief>>,
    mut contexts: bevy_egui::EguiContexts,
) -> Result {
    if contexts.ctx_mut()?.egui_wants_keyboard_input() {
        return Ok(());
    }
    if let Some(mut relief) = relief
        && keys.just_pressed(KeyCode::KeyH)
    {
        relief.raised = !relief.raised;
    }
    Ok(())
}

/// Ease the surface towards its target, and rebuild the mesh from it.
pub fn animate(time: Res<Time>, relief: Option<ResMut<Relief>>, mut meshes: ResMut<Assets<Mesh>>) {
    let Some(mut relief) = relief else { return };
    let relief = &mut *relief;
    let k = 1.0 - (-EASE * time.delta_secs()).exp();
    let lift_goal = if relief.raised { 1.0 } else { 0.0 };
    let mut moving = (relief.lift - lift_goal).abs() > 1e-3;
    relief.lift += (lift_goal - relief.lift) * k;
    for (now, goal) in relief.current.iter_mut().zip(&relief.target) {
        let delta = goal - *now;
        if delta.abs() > 1e-3 {
            *now += delta * k;
            moving = true;
        }
    }
    if !moving {
        return;
    }

    let (w, h) = (relief.width, relief.height);
    let lift = relief.lift;
    let height = |x: usize, y: usize| relief.current[y * w + x] * lift;
    let mut positions = Vec::with_capacity(w * h);
    let mut normals = Vec::with_capacity(w * h);
    for y in 0..h {
        for x in 0..w {
            positions.push(at(x, y, w, h, height(x, y)));
            // Central differences, one-sided at the edge of the box.
            let (x0, x1) = (x.saturating_sub(1), (x + 1).min(w - 1));
            let (y0, y1) = (y.saturating_sub(1), (y + 1).min(h - 1));
            let dx = (height(x1, y) - height(x0, y)) / (x1 - x0) as f32;
            let dz = (height(x, y1) - height(x, y0)) / (y1 - y0) as f32;
            normals.push(Vec3::new(-dx, 1.0, -dz).normalize().to_array());
        }
    }
    if let Some(mut mesh) = meshes.get_mut(&relief.mesh) {
        mesh.insert_attribute(Mesh::ATTRIBUTE_POSITION, positions);
        mesh.insert_attribute(Mesh::ATTRIBUTE_NORMAL, normals);
    }
}
