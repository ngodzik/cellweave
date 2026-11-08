//! One cell opened: click a cell, and a window beside the tissue shows its
//! network as it stands, every protein at its level and every link as hard
//! as it pushes, what the cell reads from where it sits, and what it does
//! with it.
//!
//! Nothing in the window is drawn for the picture's sake: every number is
//! read from the simulation, between two batches, as the census is.

use bevy::prelude::*;
use bevy::window::PrimaryWindow;
use bevy_egui::{EguiContexts, egui};
use cellweave_core::{CellId, CellInputs, Roles, SignalingNetwork};
use cellweave_engine::Fate;
use cellweave_engine::cpm::OBSTACLE;
use cellweave_engine::simulation::{DEATH_THRESHOLD, DIVISION_RATIO, QUIESCENCE_THRESHOLD};
use egui_plot::{Legend, Line, Plot, PlotPoints};

use crate::diagram::{self, Clicked, Layout};
use crate::orbit::Orbit;
use crate::paint::{self, View};
use crate::relief::Relief;
use crate::runner::Runner;

/// Pixels the pointer may move between press and release for the two to
/// be a click on a cell rather than a turn of the camera.
const CLICK: f32 = 5.0;

/// Readings kept at most; past this every other one is dropped, as the
/// panel's history is.
const HISTORY_LIMIT: usize = 1500;

/// How quickly the camera catches up with the cell it follows, per second.
const FOLLOW: f32 = 4.0;

/// How high above the tissue the halo around the open cell floats, in
/// pixels: over the tallest dome a living cell makes.
const HALO_HEIGHT: f32 = 6.5;

/// How bright the halo is, in HDR units: enough to bloom.
const HALO_GLOW: f32 = 6.0;

/// How wide the window opens: room at full size for a network of nine
/// columns, the fuller cell's.
const WIDTH: f32 = 860.0;

/// One term of the cell's network, as it stood when last read.
#[derive(Debug, Clone)]
pub struct Link {
    /// The node it moves, by index.
    pub node: usize,
    /// Its kind, as a network file names it.
    pub kind: &'static str,
    /// What it reads, by name, if anything.
    pub source: Option<String>,
    /// Its rate.
    pub rate: f64,
    /// What it adds to the node's rate of change.
    pub now: f64,
}

/// The cell as it was last read.
#[derive(Debug, Clone)]
pub struct Reading {
    /// The step it was read at.
    pub mcs: u64,
    /// Its proteins, by name, and their levels.
    pub names: Vec<String>,
    pub levels: Vec<f64>,
    /// Every term, as it pushes.
    pub links: Vec<Link>,
    /// What it reads from where it sits.
    pub inputs: CellInputs,
    /// What it is doing; `None` once it is gone.
    pub fate: Option<Fate>,
    /// Its size, what it is trying to be, and the size its lineage divides
    /// against.
    pub volume: u32,
    pub target: u32,
    pub reference: u32,
    /// Its clone.
    pub clone: Option<u32>,
    /// Its centre, in pixels.
    pub centre: Option<Vec2>,
}

/// Which cell is open, and what was read of it.
#[derive(Resource, Default)]
pub struct Inspected {
    /// The run it belongs to, by its place on screen, and the cell.
    pub cell: Option<(usize, CellId)>,
    /// A pixel clicked, waiting for its run to be back from the pool to say
    /// which cell holds it.
    clicked: Option<(usize, u32, u32)>,
    /// Keep the camera on the cell.
    pub follow: bool,
    /// The latest reading, kept once the cell is gone.
    pub reading: Option<Reading>,
    /// Every reading's levels, oldest first.
    pub history: Vec<(u64, Vec<f64>)>,
    /// The generation last read.
    seen: Option<u64>,
    /// Where a press of the left button began.
    pressed_at: Option<Vec2>,
    /// A pixel of the first run to open once it reaches a step, asked for
    /// on the command line.
    opening: Option<(u32, u32, u64)>,
    /// Where the window's left edge is, as a share of the screen's width,
    /// so that a cell followed can be kept in what it leaves uncovered.
    covered_from: Option<f32>,
}

impl Inspected {
    /// The cell open in the run at `index`, if there is one.
    pub fn in_run(&self, index: usize) -> Option<CellId> {
        self.cell.filter(|(i, _)| *i == index).map(|(_, id)| id)
    }

    /// Open the cell under pixel `(x, y)` of the first run once it reaches
    /// `step`, following it with the camera if `follow`.
    pub fn opening(self, x: u32, y: u32, step: u64, follow: bool) -> Self {
        Self {
            opening: Some((x, y, step)),
            follow,
            ..self
        }
    }

    fn close(&mut self) {
        self.cell = None;
        self.reading = None;
        self.history.clear();
        self.seen = None;
        self.follow = false;
    }
}

/// A click that does not drag picks the pixel under the pointer, on
/// whichever run it falls; the cell holding it is found once that run is
/// back from the pool.
#[allow(clippy::too_many_arguments)]
pub fn pick(
    buttons: Res<ButtonInput<MouseButton>>,
    windows: Query<&Window, With<PrimaryWindow>>,
    cameras: Query<(&Camera, &GlobalTransform), With<Camera3d>>,
    sides: Query<(&Runner, &Relief, &GlobalTransform)>,
    mut contexts: EguiContexts,
    mut inspected: ResMut<Inspected>,
) -> Result {
    let Ok(window) = windows.single() else {
        return Ok(());
    };
    let Some(cursor) = window.cursor_position() else {
        return Ok(());
    };
    let ctx = contexts.ctx_mut()?;
    if ctx.is_pointer_over_egui() || ctx.egui_wants_pointer_input() {
        inspected.pressed_at = None;
        return Ok(());
    }
    if buttons.just_pressed(MouseButton::Left) {
        inspected.pressed_at = Some(cursor);
    }
    if !buttons.just_released(MouseButton::Left) {
        return Ok(());
    }
    let Some(start) = inspected.pressed_at.take() else {
        return Ok(());
    };
    if start.distance(cursor) > CLICK {
        return Ok(());
    }
    let Ok((camera, eye)) = cameras.single() else {
        return Ok(());
    };
    let Ok(ray) = camera.viewport_to_world(eye, cursor) else {
        return Ok(());
    };
    match hit(ray, &sides) {
        Some(pixel) => inspected.clicked = Some(pixel),
        None => inspected.close(),
    }
    Ok(())
}

/// The first pixel of any run's surface the ray meets, and the run's
/// place: a march down the ray from above the tallest dome to the plane,
/// stopping where it passes under the surface.
fn hit(
    ray: Ray3d,
    sides: &Query<(&Runner, &Relief, &GlobalTransform)>,
) -> Option<(usize, u32, u32)> {
    let (origin, direction) = (ray.origin, *ray.direction);
    if direction.y >= 0.0 {
        return None;
    }
    let top = 10.0;
    let t0 = ((top - origin.y) / direction.y).max(0.0);
    let t1 = -origin.y / direction.y;
    let steps = 400;
    for k in 0..=steps {
        let point = origin + direction * (t0 + (t1 - t0) * k as f32 / steps as f32);
        for (runner, relief, place) in sides {
            let (w, h) = relief.size();
            let local = point - place.translation();
            let x = local.x + w as f32 / 2.0;
            let y = local.z + h as f32 / 2.0;
            if x < 0.0 || y < 0.0 || x >= w as f32 || y >= h as f32 {
                continue;
            }
            let (x, y) = (x as usize, y as usize);
            if point.y <= relief.height(x, y) + 0.05 {
                return Some((runner.index, x as u32, y as u32));
            }
        }
    }
    None
}

/// Read the open cell, each time its run is back from the pool with a new
/// state, and find the cell a click asked for.
pub fn read(runners: Query<&Runner>, mut inspected: ResMut<Inspected>) {
    if let Some((x, y, step)) = inspected.opening
        && runners
            .iter()
            .find(|r| r.index == 0)
            .and_then(Runner::latest)
            .is_some_and(|m| m.mcs >= step)
    {
        inspected.opening = None;
        inspected.clicked = Some((0, x, y));
    }
    if let Some((index, x, y)) = inspected.clicked {
        let Some(runner) = runners.iter().find(|r| r.index == index) else {
            inspected.clicked = None;
            return;
        };
        let Some(sim) = runner.idle() else {
            return;
        };
        inspected.clicked = None;
        let occupant = sim.lattice().occupant(x, y);
        let id = CellId(occupant);
        if occupant == 0 || occupant == OBSTACLE || sim.network(id).is_none() {
            inspected.close();
            return;
        }
        let follow = inspected.follow;
        inspected.close();
        inspected.cell = Some((index, id));
        inspected.follow = follow;
    }

    let Some((index, id)) = inspected.cell else {
        return;
    };
    let Some(runner) = runners.iter().find(|r| r.index == index) else {
        return;
    };
    if inspected
        .seen
        .is_some_and(|seen| runner.restarted_since(seen))
    {
        inspected.close();
        return;
    }
    let Some(sim) = runner.idle() else {
        return;
    };
    if inspected.seen == Some(runner.generation()) {
        return;
    }
    inspected.seen = Some(runner.generation());
    let mcs = runner.latest().map_or(0, |m| m.mcs);

    let (Some(net), Some(inputs)) = (sim.network(id), sim.inputs(id)) else {
        // Gone: the last reading stays, marked as such.
        if let Some(reading) = inspected.reading.as_mut() {
            reading.fate = None;
        }
        return;
    };
    let lattice = sim.lattice();
    let (w, _) = lattice.dims();
    let pixels = lattice.pixels_of(&[id]).pop().unwrap_or_default();
    let centre = (!pixels.is_empty()).then(|| {
        let sum = pixels.iter().fold(Vec2::ZERO, |sum, &i| {
            sum + Vec2::new((i % w as usize) as f32, (i / w as usize) as f32)
        });
        sum / pixels.len() as f32
    });
    let levels: Vec<f64> = net.state().values.iter().map(|c| c.0).collect();
    let reading = Reading {
        mcs,
        names: net.node_names().to_vec(),
        levels: levels.clone(),
        links: net
            .influences(&inputs)
            .into_iter()
            .map(|i| Link {
                node: i.node,
                kind: i.kind,
                source: i.source.map(str::to_owned),
                rate: i.rate,
                now: i.now,
            })
            .collect(),
        inputs,
        fate: sim.fate(id),
        volume: lattice.volume(id).map_or(0, |v| v.0),
        target: lattice.target_volume(id).map_or(0, |v| v.0),
        reference: net.reference_volume(),
        clone: sim.clone_of(id).map(|c| c.0),
        centre,
    };
    inspected.reading = Some(reading);
    inspected.history.push((mcs, levels));
    if inspected.history.len() > HISTORY_LIMIT {
        let last = inspected.history.len() - 1;
        let mut k = 0;
        inspected.history.retain(|_| {
            let keep = k % 2 == 0 || k == last;
            k += 1;
            keep
        });
    }
}

/// A ring of light over the open cell, so it can be found in the tissue at
/// any zoom: the white membrane alone is lost on a dome seen from afar.
#[derive(Component)]
pub struct Halo;

/// Make the halo, hidden until a cell is opened.
pub fn spawn_halo(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
) {
    commands.spawn((
        Halo,
        Mesh3d(meshes.add(Annulus::new(0.82, 1.0))),
        // Unlit, its colour is what it gives off: above one, it blooms.
        MeshMaterial3d(materials.add(StandardMaterial {
            base_color: Color::linear_rgb(HALO_GLOW, HALO_GLOW * 0.85, HALO_GLOW * 0.55),
            unlit: true,
            ..default()
        })),
        // An annulus lies in the plane facing +Z; the tissue faces +Y.
        Transform::from_rotation(Quat::from_rotation_x(-std::f32::consts::FRAC_PI_2)),
        Visibility::Hidden,
    ));
}

/// Put the halo over the open cell, as wide as the cell, easing towards
/// where the cell has moved since the last batch.
pub fn halo(
    time: Res<Time>,
    inspected: Res<Inspected>,
    sides: Query<(&Runner, &Relief, &GlobalTransform)>,
    mut halos: Query<(&mut Transform, &mut Visibility), With<Halo>>,
) {
    let place = inspected.cell.and_then(|(index, _)| {
        let reading = inspected.reading.as_ref()?;
        reading.fate?;
        let centre = reading.centre?;
        let (_, relief, side) = sides.iter().find(|(r, _, _)| r.index == index)?;
        let (w, h) = relief.size();
        let at = side.translation()
            + Vec3::new(
                centre.x + 0.5 - w as f32 / 2.0,
                HALO_HEIGHT,
                centre.y + 0.5 - h as f32 / 2.0,
            );
        let radius = (reading.volume as f32 / std::f32::consts::PI).sqrt() + 1.0;
        Some((at, radius))
    });
    let k = 1.0 - (-2.0 * FOLLOW * time.delta_secs()).exp();
    for (mut transform, mut visibility) in &mut halos {
        let Some((at, radius)) = place else {
            *visibility = Visibility::Hidden;
            continue;
        };
        if *visibility == Visibility::Hidden {
            // Appearing, it starts where the cell is rather than flying in.
            transform.translation = at;
            transform.scale = Vec3::splat(radius);
        }
        *visibility = Visibility::Visible;
        let (now, size) = (transform.translation, transform.scale);
        transform.translation = now.lerp(at, k);
        transform.scale = size.lerp(Vec3::splat(radius), k);
    }
}

/// Keep the camera on the cell while asked to, and ease it back to the
/// centre of the tissue once not.
///
/// The cell is kept in the middle of what the window leaves of the screen
/// rather than in the middle of the screen, where the window would hide it:
/// the camera looks that much to the right of it.
pub fn follow(
    time: Res<Time>,
    inspected: Res<Inspected>,
    sides: Query<(&Runner, &Relief, &GlobalTransform)>,
    mut orbit: ResMut<Orbit>,
    mut camera: Query<(&mut Transform, &Projection), With<Camera3d>>,
) {
    let aside = camera
        .iter()
        .next()
        .and_then(|(transform, projection)| {
            let Projection::Perspective(lens) = projection else {
                return None;
            };
            let from = inspected.covered_from?.clamp(0.0, 1.0);
            let half_width = orbit.distance() * (lens.fov / 2.0).tan() * lens.aspect_ratio;
            Some(transform.right() * (1.0 - from) * half_width)
        })
        .unwrap_or(Vec3::ZERO);
    let goal = inspected
        .follow
        .then(|| {
            let (index, _) = inspected.cell?;
            let centre = inspected.reading.as_ref()?.centre?;
            let (_, relief, place) = sides.iter().find(|(r, _, _)| r.index == index)?;
            let (w, h) = relief.size();
            Some(
                place.translation()
                    + Vec3::new(
                        centre.x + 0.5 - w as f32 / 2.0,
                        0.0,
                        centre.y + 0.5 - h as f32 / 2.0,
                    ),
            )
        })
        .flatten()
        .map_or(Vec3::ZERO, |cell| cell + aside);
    let k = 1.0 - (-FOLLOW * time.delta_secs()).exp();
    if orbit.ease_focus(goal, k) {
        for (mut transform, _) in &mut camera {
            *transform = orbit.transform();
        }
    }
}

/// The window, when a cell is open.
pub fn window(
    mut contexts: EguiContexts,
    mut inspected: ResMut<Inspected>,
    runners: Query<&Runner>,
    mut view: ResMut<View>,
    mut orbit: ResMut<Orbit>,
) -> Result {
    let Some((index, id)) = inspected.cell else {
        return Ok(());
    };
    let Some(runner) = runners.iter().find(|r| r.index == index) else {
        return Ok(());
    };
    let roles = runner.experiment.spec.roles.clone();
    let ctx = contexts.ctx_mut()?;
    let mut open = true;
    let mut chosen = match *view {
        View::Protein(i) => Some(i),
        _ => None,
    };
    let mut follow = inspected.follow;
    let inspected_ref = &*inspected;
    let screen = ctx.content_rect();
    let right = screen.right();
    let shown = egui::Window::new(format!("Inside cell {}", id.0))
        .id(egui::Id::new("inside a cell"))
        .default_pos([right - WIDTH - 24.0, 12.0])
        .default_width(WIDTH)
        .resizable(false)
        .open(&mut open)
        .show(ctx, |ui| {
            let Some(reading) = inspected_ref.reading.as_ref() else {
                ui.label("Reading the cell…");
                return;
            };
            header(ui, reading, runner);
            ui.checkbox(&mut follow, "follow it with the camera");
            ui.separator();
            let badges = badges(reading, roles.as_ref());
            let layout = Layout::of(reading.names.len(), &reading.names, &reading.links);
            let reads = |name: &str| read_input(&reading.inputs, name);
            let explain = |n: usize| explain(reading, n, roles.as_ref());
            if let Clicked::Node(n) = diagram::draw(
                ui,
                &layout,
                &reading.names,
                &reading.levels,
                &reading.links,
                &reads,
                &badges,
                chosen,
                &explain,
            ) {
                chosen = if chosen == Some(n) { None } else { Some(n) };
            }
            ui.horizontal_wrapped(|ui| {
                legend(ui, egui::Color32::from_rgb(90, 205, 150), "switches on");
                legend(ui, egui::Color32::from_rgb(240, 95, 90), "switches off");
                legend(ui, egui::Color32::from_rgb(190, 60, 120), "damages");
                ui.weak("· dashed: by its absence · thicker: pushing harder");
            });
            ui.separator();
            decisions(ui, reading, roles.as_ref());
            ui.separator();
            history(ui, &inspected_ref.history, reading, roles.as_ref(), chosen);
            ui.separator();
            ui.weak(
                "Click a protein to colour the tissue by it. Levels run from 0 to 1, as the model computes them. \
                 The links come from the literature; their rates are chosen, not measured.",
            );
        });
    inspected.covered_from = shown.map(|s| s.response.rect.left() / screen.width().max(1.0));
    if !open {
        inspected.close();
        return Ok(());
    }
    if follow && !inspected.follow {
        // Coming closer once, when asked to follow, so the cell can be seen.
        *orbit = orbit.closer_than(0.35);
    }
    inspected.follow = follow;
    let wanted = match chosen {
        Some(i) => View::Protein(i),
        None if matches!(*view, View::Protein(_)) => View::Fate,
        None => *view,
    };
    if wanted != *view {
        *view = wanted;
    }
    Ok(())
}

/// Which cell, what it is doing, and how big it is against what it wants.
fn header(ui: &mut egui::Ui, reading: &Reading, runner: &Runner) {
    let (fate, colour) = match reading.fate {
        Some(Fate::Cycling) => ("cycling: growing to divide", paint::CYCLING),
        Some(Fate::Quiescent) => ("arrested: alive, not dividing", paint::QUIESCENT),
        Some(Fate::Necrotic) => ("necrotic: dead, giving up its pixels", paint::NECROTIC),
        Some(Fate::Differentiating) => {
            ("differentiating: leaving the layer", paint::DIFFERENTIATING)
        }
        Some(Fate::Apoptotic) => ("apoptotic: breaking up", paint::APOPTOTIC),
        None => ("gone: its last reading", paint::NECROTIC),
    };
    ui.horizontal(|ui| {
        let (rect, _) = ui.allocate_exact_size(egui::vec2(12.0, 12.0), egui::Sense::hover());
        ui.painter().rect_filled(rect, 2.0, swatch(colour));
        ui.strong(fate);
    });
    let divides = (DIVISION_RATIO * f64::from(reading.reference)).round();
    ui.label(format!(
        "{} pixels, asking for {}; divides at {divides:.0}, if cycling · clone {} · {} · step {}",
        reading.volume,
        reading.target,
        reading.clone.map_or("none".to_owned(), |c| c.to_string()),
        runner.name,
        reading.mcs,
    ));
}

/// A short label under each protein that answers one of the questions the
/// simulation asks a cell, lit while that answer is the one it gives.
fn badges(reading: &Reading, roles: Option<&Roles>) -> Vec<Option<(String, bool)>> {
    let mut out = vec![None; reading.names.len()];
    let Some(roles) = roles else {
        return out;
    };
    let level = |n: usize| reading.levels.get(n).copied().unwrap_or(0.0);
    let mut set = |name: &str, badge: &str, on: &dyn Fn(f64) -> bool| {
        if let Some(n) = reading.names.iter().position(|x| x == name) {
            let previous: Option<(String, bool)> = out[n].take();
            let text = match previous {
                Some((before, _)) => format!("{before} · {badge}"),
                None => badge.to_owned(),
            };
            out[n] = Some((text, on(level(n))));
        }
    };
    set(&roles.growth, "size", &|_| true);
    set(&roles.survival, "alive", &|v| v >= DEATH_THRESHOLD);
    set(&roles.hypoxia_response, "arrest", &|v| {
        v >= QUIESCENCE_THRESHOLD
    });
    if let Some(name) = &roles.apoptosis {
        set(name, "apoptosis", &|v| v > 0.01);
    }
    if let Some(name) = &roles.angiogenic {
        set(name, "vessels", &|v| v > 0.05);
    }
    if let Some(name) = &roles.differentiation {
        set(name, "leaves", &|v| v > 0.01);
    }
    if let Some(name) = &roles.evasion {
        set(name, "hides", &|v| v > 0.05);
    }
    out
}

/// What the cell does with its levels, each answer with the protein and
/// the threshold it was read against.
fn decisions(ui: &mut egui::Ui, reading: &Reading, roles: Option<&Roles>) {
    let Some(roles) = roles else {
        return;
    };
    let level = |name: &str| {
        reading
            .names
            .iter()
            .position(|n| n == name)
            .and_then(|n| reading.levels.get(n).copied())
            .unwrap_or(0.0)
    };
    let name = |n: &str| diagram::label(n);
    egui::Grid::new("decisions").num_columns(2).show(ui, |ui| {
        let arrest = level(&roles.hypoxia_response);
        ui.label("divide?");
        ui.label(if arrest < QUIESCENCE_THRESHOLD {
            format!(
                "yes: {} at {arrest:.3}, below the {QUIESCENCE_THRESHOLD} that stops the cycle",
                name(&roles.hypoxia_response)
            )
        } else {
            format!(
                "no: {} at {arrest:.3}, at or above the {QUIESCENCE_THRESHOLD} that stops the cycle",
                name(&roles.hypoxia_response)
            )
        });
        ui.end_row();
        let growth = level(&roles.growth);
        ui.label("how big?");
        ui.label(format!(
            "{} at {growth:.2} asks for {} pixels, against a reference of {}",
            name(&roles.growth),
            reading.target,
            reading.reference
        ));
        ui.end_row();
        let survival = level(&roles.survival);
        ui.label("alive?");
        ui.label(if survival >= DEATH_THRESHOLD {
            format!(
                "yes: {} at {survival:.3}, above the {DEATH_THRESHOLD} below which it dies",
                name(&roles.survival)
            )
        } else {
            format!(
                "no: {} at {survival:.3}, fell below {DEATH_THRESHOLD}",
                name(&roles.survival)
            )
        });
        ui.end_row();
        if let Some(apoptosis) = &roles.apoptosis {
            let chance = level(apoptosis);
            ui.label("self-destruct?");
            ui.label(format!(
                "{} at {chance:.3}: a {:.1}% chance each step",
                name(apoptosis),
                100.0 * chance
            ));
            ui.end_row();
        }
        ui.label("reads");
        let mut reads = Vec::new();
        for input in cellweave_core::network::INPUT_NAMES {
            if reading
                .links
                .iter()
                .any(|l| l.source.as_deref() == Some(input))
            {
                reads.push(format!(
                    "{} {:.2}",
                    diagram::label(input),
                    read_input(&reading.inputs, input)
                ));
            }
        }
        ui.label(reads.join(" · "));
        ui.end_row();
    });
}

/// The levels over the cell's life since it was opened: the proteins that
/// decide, and the one the tissue is coloured by.
fn history(
    ui: &mut egui::Ui,
    history: &[(u64, Vec<f64>)],
    reading: &Reading,
    roles: Option<&Roles>,
    chosen: Option<usize>,
) {
    let mut shown: Vec<usize> = Vec::new();
    if let Some(roles) = roles {
        for name in [
            Some(&roles.growth),
            Some(&roles.hypoxia_response),
            Some(&roles.survival),
            roles.apoptosis.as_ref(),
        ]
        .into_iter()
        .flatten()
        {
            if let Some(n) = reading.names.iter().position(|x| x == name)
                && !shown.contains(&n)
            {
                shown.push(n);
            }
        }
    }
    if let Some(n) = chosen
        && !shown.contains(&n)
    {
        shown.push(n);
    }
    const COLOURS: [egui::Color32; 6] = [
        egui::Color32::from_rgb(90, 205, 150),
        egui::Color32::from_rgb(255, 170, 60),
        egui::Color32::from_rgb(120, 170, 255),
        egui::Color32::from_rgb(200, 120, 240),
        egui::Color32::WHITE,
        egui::Color32::from_rgb(240, 95, 90),
    ];
    Plot::new("inside a cell, over time")
        .height(110.0)
        .legend(Legend::default())
        .allow_drag(false)
        .allow_zoom(false)
        .allow_scroll(false)
        .allow_boxed_zoom(false)
        .include_y(0.0)
        .include_y(1.0)
        .show(ui, |plot| {
            for (k, &n) in shown.iter().enumerate() {
                let points: PlotPoints = history
                    .iter()
                    .filter_map(|(mcs, levels)| Some([*mcs as f64, *levels.get(n)?]))
                    .collect();
                let colour = if Some(n) == chosen {
                    egui::Color32::WHITE
                } else {
                    COLOURS[k % COLOURS.len()]
                };
                plot.line(Line::new(diagram::label(&reading.names[n]), points).color(colour));
            }
        });
}

/// What hovering over a protein says: what it is for, if it answers one of
/// the simulation's questions, and every term that moves it, as it pushes.
fn explain(reading: &Reading, n: usize, roles: Option<&Roles>) -> String {
    let name = &reading.names[n];
    let mut lines = vec![format!(
        "{}: {:.3}",
        diagram::label(name),
        reading.levels.get(n).copied().unwrap_or(0.0)
    )];
    if let Some(roles) = roles {
        let role = [
            (
                Some(&roles.growth),
                "its level sets the size the cell asks for",
            ),
            (
                Some(&roles.survival),
                "the cell dies when it falls below 0.25",
            ),
            (
                Some(&roles.hypoxia_response),
                "the cell stops dividing when it reaches 0.45",
            ),
            (
                roles.apoptosis.as_ref(),
                "its level is the chance each step that the cell starts apoptosis",
            ),
            (
                roles.angiogenic.as_ref(),
                "its level is what the cell secretes as a call for vessels",
            ),
            (
                roles.differentiation.as_ref(),
                "its level is the chance each step that the cell leaves the layer",
            ),
            (
                roles.evasion.as_ref(),
                "its level shields the cell from immune cells",
            ),
        ];
        for (who, what) in role {
            if who == Some(name) {
                lines.push(format!("Role: {what}."));
            }
        }
    }
    for link in reading.links.iter().filter(|l| l.node == n) {
        let what = match (link.kind, link.source.as_deref()) {
            ("activated_by", Some(s)) => format!("switched on by {}", diagram::label(s)),
            ("inhibited_by", Some(s)) => format!("switched off by {}", diagram::label(s)),
            ("activated_by_absence_of", Some(s)) => {
                format!("switched on by the absence of {}", diagram::label(s))
            }
            ("produced_by", Some(s)) => format!("made in proportion to {}", diagram::label(s)),
            ("damaged_below", Some(s)) => format!("damaged when {} runs short", diagram::label(s)),
            ("damaged_above", Some(s)) => format!("damaged when {} runs high", diagram::label(s)),
            ("decays", _) => "turned over".to_owned(),
            ("relaxes_to", _) => "held near a set point".to_owned(),
            (kind, _) => kind.to_owned(),
        };
        lines.push(format!(
            "{what}, rate {:.2}: {:+.3} now",
            link.rate, link.now
        ));
    }
    lines.join("\n")
}

/// The value of an input, by the name a network file gives it.
fn read_input(inputs: &CellInputs, name: &str) -> f64 {
    match name {
        "oxygen" => inputs.o2.0,
        "egf" => inputs.egf.0,
        "crowding" => inputs.crowding.0,
        "stress" => inputs.stress.0,
        "waste" => inputs.waste.0,
        "necrotic_factor" => inputs.necrotic_factor.0,
        _ => 0.0,
    }
}

fn legend(ui: &mut egui::Ui, colour: egui::Color32, what: &str) {
    let (rect, _) = ui.allocate_exact_size(egui::vec2(16.0, 10.0), egui::Sense::hover());
    ui.painter().line_segment(
        [rect.left_center(), rect.right_center()],
        egui::Stroke::new(3.0, colour),
    );
    ui.label(what);
}

fn swatch(c: LinearRgba) -> egui::Color32 {
    let [r, g, b, _] = paint::bytes(c);
    egui::Color32::from_rgb(r, g, b)
}
