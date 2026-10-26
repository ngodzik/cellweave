//! The light the tissue gives off, separate from its colour under the lamp.
//!
//! An image the size of the canvas, used as the surface's emission: living
//! cells glow from within, dark at their membranes so the light does not
//! wash the outlines out, the vessel burns, and the events of a step flash
//! and fade. A cell born in the last
//! batch flashes white; a cell that has just died flares red and goes out.
//! The simulation moves in jumps of a few steps; the glow is redrawn every
//! frame, so what happened in a jump is still visible a second later.

use bevy::asset::RenderAssetUsages;
use bevy::prelude::*;
use bevy::render::render_resource::{Extent3d, TextureDimension, TextureFormat};
use cellweave_core::{CellId, DiffusionField};
use cellweave_engine::cpm::OBSTACLE;
use cellweave_engine::simulation::Running;
use cellweave_engine::{Fate, Simulation};

use crate::paint::{self, View};
use crate::runner::Runner;

/// How much of a cell's colour it gives off as light, by fate. A cycling
/// cell is the most alive; an arrested one smoulders; a dead one is dark.
const CYCLING_GLOW: f32 = 0.06;
const QUIESCENT_GLOW: f32 = 0.04;
const VESSEL_GLOW: f32 = 0.35;

/// What the medium gives off where it is full of oxygen: a faint blue haze
/// that shows where the supply reaches.
const MEDIUM_GLOW: LinearRgba = LinearRgba::rgb(0.0, 0.004, 0.012);

/// A newborn's flash, and how long it takes to fade, in seconds.
const BIRTH: LinearRgba = LinearRgba::rgb(0.30, 0.45, 0.42);
const BIRTH_FADE: f32 = 0.22;
/// A death's flare, and how long it takes to go out, in seconds.
const DEATH: LinearRgba = LinearRgba::rgb(0.45, 0.07, 0.015);
const DEATH_FADE: f32 = 0.35;

/// The emission image and what it is drawn from.
#[derive(Resource)]
pub struct Glow {
    /// Written every frame something is fading.
    pub image: Handle<Image>,
    /// The material showing it, touched so the new image is picked up.
    pub material: Handle<StandardMaterial>,
    width: u32,
    /// Who holds each pixel, as of the last state drawn.
    occupant: Vec<u32>,
    /// Which sides of each pixel are a membrane.
    edges: Vec<[bool; 4]>,
    /// What each pixel gives off before any flash: the cell's steady glow,
    /// or the medium's haze.
    steady: Vec<LinearRgba>,
    /// Per cell identifier: how bright its birth flash still is, from one.
    birth: Vec<f32>,
    /// Per cell identifier: how bright its death flare still is, from one.
    death: Vec<f32>,
    /// Per cell identifier: its fate in the last state drawn.
    fates: Vec<Option<Fate>>,
    /// The generation last drawn, and the view it was drawn in.
    seen: Option<(u64, View)>,
    /// Something is still fading, so the image changes every frame.
    fading: bool,
}

impl Glow {
    /// A dark glow for a lattice this size.
    pub fn new(
        image: Handle<Image>,
        material: Handle<StandardMaterial>,
        width: u32,
        height: u32,
    ) -> Self {
        let n = (width * height) as usize;
        Self {
            image,
            material,
            width,
            occupant: vec![0; n],
            edges: vec![[false; 4]; n],
            steady: vec![LinearRgba::BLACK; n],
            birth: Vec::new(),
            death: Vec::new(),
            fates: Vec::new(),
            seen: None,
            fading: false,
        }
    }
}

/// A black emission image, the size of the canvas.
pub fn image(width: u32, height: u32) -> Image {
    Image::new_fill(
        Extent3d {
            width: width * paint::SCALE,
            height: height * paint::SCALE,
            depth_or_array_layers: 1,
        },
        TextureDimension::D2,
        &[0, 0, 0, 255],
        TextureFormat::Rgba8UnormSrgb,
        RenderAssetUsages::MAIN_WORLD | RenderAssetUsages::RENDER_WORLD,
    )
}

/// Take in a new state: what each pixel gives off, which cells were just
/// born, and which just died.
pub fn retarget(runner: Res<Runner>, view: Res<View>, glow: Option<ResMut<Glow>>) {
    let Some(mut glow) = glow else { return };
    let Some(sim) = runner.idle() else { return };
    let now = (runner.generation(), *view);
    if glow.seen == Some(now) {
        return;
    }
    let restarted = glow
        .seen
        .is_none_or(|(seen, _)| runner.restarted_since(seen));
    glow.seen = Some(now);
    let glow = &mut *glow;

    let lattice = sim.lattice();
    let cells = lattice.cell_count() + 1;
    if restarted {
        glow.birth.clear();
        glow.death.clear();
        glow.fates.clear();
    }
    let fates: Vec<Option<Fate>> = (0..cells as u32).map(|id| sim.fate(CellId(id))).collect();
    // A first state, or a restart, has no events: every cell in it would
    // otherwise flash as newborn.
    let known = if glow.fates.is_empty() {
        cells
    } else {
        glow.fates.len()
    };
    glow.birth.resize(cells, 0.0);
    glow.death.resize(cells, 0.0);
    for id in known..cells {
        glow.birth[id] = 1.0;
    }
    for (id, fate) in fates.iter().enumerate() {
        let was_alive = matches!(
            glow.fates.get(id),
            Some(Some(Fate::Cycling | Fate::Quiescent))
        );
        if was_alive && *fate == Some(Fate::Necrotic) {
            glow.death[id] = 1.0;
        }
    }
    glow.fates = fates;

    let per_cell: Vec<LinearRgba> = (0..cells as u32)
        .map(|id| steady(sim, CellId(id), *view, glow.fates[id as usize]))
        .collect();
    let oxygen = sim.oxygen();
    let (w, h) = lattice.dims();
    for y in 0..h {
        for x in 0..w {
            let at = (y * glow.width + x) as usize;
            let id = lattice.occupant(x, y);
            glow.occupant[at] = id;
            glow.edges[at] = paint::membrane(lattice, x, y);
            glow.steady[at] = match id {
                0 => MEDIUM_GLOW * oxygen.concentration_at(x, y) as f32,
                OBSTACLE => paint::VESSEL * VESSEL_GLOW,
                _ => per_cell
                    .get(id as usize)
                    .copied()
                    .unwrap_or(LinearRgba::BLACK),
            };
        }
    }
    glow.fading = true;
}

/// What a cell gives off before any flash.
fn steady(sim: &Simulation<Running>, id: CellId, view: View, fate: Option<Fate>) -> LinearRgba {
    let strength = match fate {
        Some(Fate::Cycling) => CYCLING_GLOW,
        Some(Fate::Quiescent) => QUIESCENT_GLOW,
        Some(Fate::Necrotic) | None => return LinearRgba::BLACK,
    };
    paint::cell_colour(sim, id, view) * strength
}

/// Fade every flash a frame's worth and redraw the image, while anything is
/// still changing.
pub fn animate(
    time: Res<Time>,
    glow: Option<ResMut<Glow>>,
    mut images: ResMut<Assets<Image>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    let Some(mut glow) = glow else { return };
    if !glow.fading {
        return;
    }
    let glow = &mut *glow;
    let dt = time.delta_secs();
    let (birth_k, death_k) = ((-dt / BIRTH_FADE).exp(), (-dt / DEATH_FADE).exp());
    let mut still = false;
    for b in &mut glow.birth {
        *b *= birth_k;
        if *b < 0.01 {
            *b = 0.0;
        }
        still |= *b > 0.0;
    }
    for d in &mut glow.death {
        *d *= death_k;
        if *d < 0.01 {
            *d = 0.0;
        }
        still |= *d > 0.0;
    }

    let Some(mut image) = images.get_mut(&glow.image) else {
        return;
    };
    let Some(data) = image.data.as_mut() else {
        return;
    };
    let dark = [0, 0, 0, 255];
    for (at, (&id, &base)) in glow.occupant.iter().zip(&glow.steady).enumerate() {
        let id = id as usize;
        let flash = glow.birth.get(id).copied().unwrap_or(0.0);
        let flare = glow.death.get(id).copied().unwrap_or(0.0);
        let colour = paint::bytes(base + BIRTH * flash + DEATH * flare);
        let at = at as u32;
        let xy = (at % glow.width, at / glow.width);
        paint::block(data, glow.width, xy, colour, dark, glow.edges[at as usize]);
    }
    let _ = materials.get_mut(&glow.material);
    // One last redraw after the last flash has faded, then rest.
    glow.fading = still;
}
