//! A tissue grown in a volume, drawn as the voxels it is made of: every face
//! between a cell and the medium, and, with the cut on, every face of the
//! plane through the centre of the box, which opens the tissue like a fruit
//! cut in half and shows its core.
//!
//! The mesh is rebuilt from the lattice whenever a batch comes back from the
//! pool, or the view or the cut changes. Faces inside the tissue are never
//! drawn, so the mesh holds the surface and the cut, not the bulk.

use bevy::asset::RenderAssetUsages;
use bevy::mesh::{Indices, PrimitiveTopology};
use bevy::prelude::*;
use bevy_egui::EguiContexts;
use cellweave_core::CellId;
use cellweave_engine::Simulation;
use cellweave_engine::cpm::OBSTACLE;
use cellweave_engine::simulation::Running;

use crate::paint::{self, View};
use crate::runner::Runner;

/// How much darker a voxel on a cell's boundary is than the cell's inside,
/// so that the cells of a cut read as cells and not as one block. Milder
/// than a membrane in the plane, since a voxel is a seventh of a cell
/// across where a texel is a sixtieth.
const BOUNDARY: f32 = 0.55;

/// One run's tissue in a volume.
#[derive(Component)]
pub struct Volume {
    /// The mesh rebuilt from the lattice.
    pub mesh: Handle<Mesh>,
    /// Whether the half of the box nearest the camera's first position is
    /// cut away.
    pub cut: bool,
    /// The generation, view and cut the mesh was last built from.
    seen: Option<(u64, View, bool)>,
}

impl Volume {
    /// An empty volume, built once the first state arrives, cut open.
    pub fn new(mesh: Handle<Mesh>) -> Self {
        Self {
            mesh,
            cut: true,
            seen: None,
        }
    }
}

/// An empty mesh, kept on the CPU as well as the GPU since it is rebuilt.
pub fn empty() -> Mesh {
    Mesh::new(
        PrimitiveTopology::TriangleList,
        RenderAssetUsages::MAIN_WORLD | RenderAssetUsages::RENDER_WORLD,
    )
    .with_inserted_attribute(Mesh::ATTRIBUTE_POSITION, Vec::<[f32; 3]>::new())
    .with_inserted_attribute(Mesh::ATTRIBUTE_NORMAL, Vec::<[f32; 3]>::new())
    .with_inserted_attribute(Mesh::ATTRIBUTE_COLOR, Vec::<[f32; 4]>::new())
    .with_inserted_indices(Indices::U32(Vec::new()))
}

/// C opens and closes the cut.
pub fn toggle(
    keys: Res<ButtonInput<KeyCode>>,
    mut volumes: Query<&mut Volume>,
    mut contexts: EguiContexts,
) -> Result {
    if contexts.ctx_mut()?.egui_wants_keyboard_input() || !keys.just_pressed(KeyCode::KeyC) {
        return Ok(());
    }
    for mut volume in &mut volumes {
        volume.cut = !volume.cut;
    }
    Ok(())
}

/// Rebuild each run's mesh when its simulation moved, or the view or the
/// cut changed.
pub fn rebuild(
    mut sides: Query<(&Runner, &mut Volume)>,
    view: Res<View>,
    mut meshes: ResMut<Assets<Mesh>>,
) {
    for (runner, mut volume) in &mut sides {
        let Some(sim) = runner.idle() else { continue };
        let now = (runner.generation(), *view, volume.cut);
        if volume.seen == Some(now) {
            continue;
        }
        let Some(mut mesh) = meshes.get_mut(&volume.mesh) else {
            continue;
        };
        *mesh = build(sim, *view, volume.cut);
        volume.seen = Some(now);
    }
}

/// One of the six faces of a voxel: the step on the lattice to the voxel
/// beyond it, and the world's normal and the two edges of the face, `u`
/// across `v` along the normal so that the corners go round
/// counter-clockwise seen from outside.
struct Face {
    step: [i64; 3],
    normal: [f32; 3],
    u: [f32; 3],
    v: [f32; 3],
}

/// The six faces. The lattice's rows run along the world's Z and its slices
/// up its Y, so a tissue grown in a volume stands as it would in a dish.
const FACES: [Face; 6] = [
    Face {
        step: [1, 0, 0],
        normal: [1.0, 0.0, 0.0],
        u: [0.0, 1.0, 0.0],
        v: [0.0, 0.0, 1.0],
    },
    Face {
        step: [-1, 0, 0],
        normal: [-1.0, 0.0, 0.0],
        u: [0.0, 0.0, 1.0],
        v: [0.0, 1.0, 0.0],
    },
    Face {
        step: [0, 0, 1],
        normal: [0.0, 1.0, 0.0],
        u: [0.0, 0.0, 1.0],
        v: [1.0, 0.0, 0.0],
    },
    Face {
        step: [0, 0, -1],
        normal: [0.0, -1.0, 0.0],
        u: [1.0, 0.0, 0.0],
        v: [0.0, 0.0, 1.0],
    },
    Face {
        step: [0, 1, 0],
        normal: [0.0, 0.0, 1.0],
        u: [1.0, 0.0, 0.0],
        v: [0.0, 1.0, 0.0],
    },
    Face {
        step: [0, -1, 0],
        normal: [0.0, 0.0, -1.0],
        u: [0.0, 1.0, 0.0],
        v: [1.0, 0.0, 0.0],
    },
];

/// The surface of everything the lattice holds, and with `cut` the face of
/// the plane through the middle of its rows, coloured as the plane would be
/// in `view`.
fn build(sim: &Simulation<Running>, view: View, cut: bool) -> Mesh {
    let lattice = sim.lattice();
    let (w, h, d) = lattice.dims3();
    let (wi, hi, di) = (i64::from(w), i64::from(h), i64::from(d));
    // The half of the rows nearest the camera as it starts out is cut away.
    let kept = |y: i64| !cut || y < hi / 2;
    let occupant = |x: i64, y: i64, z: i64| -> u32 {
        if x < 0 || y < 0 || z < 0 || x >= wi || y >= hi || z >= di || !kept(y) {
            return 0;
        }
        lattice.occupant_at(((z * hi + y) * wi + x) as usize)
    };

    let cells: Vec<LinearRgba> = (0..=lattice.cell_count() as u32)
        .map(|id| paint::cell_colour(sim, CellId(id), view))
        .collect();
    let field = |i: usize| -> Option<f32> {
        match view {
            View::Oxygen => Some(sim.oxygen().value_at(i) as f32),
            View::Stress => sim.stress().map(|f| f.value_at(i) as f32),
            _ => None,
        }
    };

    let mut positions: Vec<[f32; 3]> = Vec::new();
    let mut normals: Vec<[f32; 3]> = Vec::new();
    let mut colours: Vec<[f32; 4]> = Vec::new();
    let mut indices: Vec<u32> = Vec::new();
    let centre = [wi as f32 / 2.0, di as f32 / 2.0, hi as f32 / 2.0];

    for (i, &id) in lattice.grid.iter().enumerate() {
        if id == 0 {
            continue;
        }
        let (x, y, z) = lattice.position(i);
        let (x, y, z) = (i64::from(x), i64::from(y), i64::from(z));
        if !kept(y) {
            continue;
        }
        let open: Vec<usize> = (0..FACES.len())
            .filter(|&f| {
                let [dx, dy, dz] = FACES[f].step;
                occupant(x + dx, y + dy, z + dz) == 0
            })
            .collect();
        if open.is_empty() {
            continue;
        }
        let colour = if id == OBSTACLE {
            paint::VESSEL
        } else {
            let base = match field(i) {
                Some(level) => paint::ramp(level),
                None => cells.get(id as usize).copied().unwrap_or(paint::VESSEL),
            };
            // A voxel touching another cell is on this cell's boundary.
            let boundary = FACES.iter().any(|face| {
                let [dx, dy, dz] = face.step;
                let other = occupant(x + dx, y + dy, z + dz);
                other != 0 && other != id
            });
            if boundary { base * BOUNDARY } else { base }
        };
        let rgba = colour.to_f32_array();
        // The voxel's centre in the world: rows along Z, slices up Y.
        let at = [
            x as f32 + 0.5 - centre[0],
            z as f32 + 0.5 - centre[1],
            y as f32 + 0.5 - centre[2],
        ];
        for f in open {
            let Face {
                normal: n, u, v, ..
            } = FACES[f];
            let base = positions.len() as u32;
            for (su, sv) in [(-1.0, -1.0), (1.0, -1.0), (1.0, 1.0), (-1.0, 1.0)] {
                positions.push([
                    at[0] + 0.5 * (n[0] + su * u[0] + sv * v[0]),
                    at[1] + 0.5 * (n[1] + su * u[1] + sv * v[1]),
                    at[2] + 0.5 * (n[2] + su * u[2] + sv * v[2]),
                ]);
                normals.push(n);
                colours.push(rgba);
            }
            indices.extend_from_slice(&[base, base + 1, base + 2, base, base + 2, base + 3]);
        }
    }

    Mesh::new(
        PrimitiveTopology::TriangleList,
        RenderAssetUsages::MAIN_WORLD | RenderAssetUsages::RENDER_WORLD,
    )
    .with_inserted_attribute(Mesh::ATTRIBUTE_POSITION, positions)
    .with_inserted_attribute(Mesh::ATTRIBUTE_NORMAL, normals)
    .with_inserted_attribute(Mesh::ATTRIBUTE_COLOR, colours)
    .with_inserted_indices(Indices::U32(indices))
}
