//! cellweave-viewer: watch a simulation run.
//!
//! Usage:
//!
//! ```text
//! cellweave-viewer [--config <path>] [--seed <n>] [--view fate|oxygen|clone|stress]
//!                  [--beside <path> | --undo <clone>,<clone>... [--undo-at <step>]]
//! ```
//!
//! With `--beside`, a second configuration runs to the right of the first,
//! under the first's seed. With `--undo`, the right is the same run with the
//! named clones' mutations undone at one step, by default the step its lesion
//! appears at: the two are the same run, to the last bit, until then. Either
//! way the two are stepped together, and the panel says when they part.
//!
//! It runs the same [`cellweave_engine::Simulation`] the command line runs,
//! under the same seed, so what it shows at a given step is what the command
//! line writes at that step. The tissue is drawn as a surface in relief, in
//! HDR: each cell a dome coloured by what it is doing, lit from the side so
//! its shape casts a shadow, and glowing from within while it lives. Births
//! flash and deaths flare, so the structure a run produces, a starving core
//! or a cord around its vessel, is visible at a glance and so is its motion.
//!
//! Space pauses, R restarts, H flattens the relief, A turns the camera on
//! its own, and the mouse turns and zooms it.

#![deny(unsafe_code)]

mod glow;
mod orbit;
mod paint;
mod panel;
mod relief;
mod runner;

use anyhow::{Result, anyhow};
use bevy::camera::visibility::NoFrustumCulling;
use bevy::core_pipeline::tonemapping::Tonemapping;
use bevy::light::CascadeShadowConfigBuilder;
use bevy::post_process::bloom::Bloom;
use bevy::prelude::*;
use bevy_egui::{EguiPlugin, EguiPrimaryContextPass};
use cellweave_core::{CloneId, NetworkSpec};
use cellweave_engine::Scenario;
use cellweave_engine::simulation::Reversion;
use cellweave_io::SimConfig;

/// What is being watched: the scenario and the hypothesis, kept so that a
/// restart starts the very same run again. The first run on screen is also a
/// resource, the one the panel describes.
#[derive(Resource, Clone)]
pub struct Experiment {
    /// The configuration file it came from, for the title and the panel.
    pub label: String,
    /// The box, the oxygen source, the seed.
    pub scenario: Scenario,
    /// The network every cell carries.
    pub spec: NetworkSpec,
}

fn main() -> Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let config_path = args
        .windows(2)
        .find(|w| w[0] == "--config")
        .map(|w| w[1].as_str());

    let mut config = match config_path {
        Some(p) => SimConfig::from_file(p).map_err(|e| anyhow!("loading config {p}: {e}"))?,
        None => SimConfig::default(),
    };
    // --seed, as on the command line, so a history found by a campaign can be
    // watched.
    if let Some(w) = args.windows(2).find(|w| w[0] == "--seed") {
        config.seed = w[1]
            .parse()
            .map_err(|e| anyhow!("--seed must be an integer: {e}"))?;
    }
    let view = match args.windows(2).find(|w| w[0] == "--view") {
        None => paint::View::default(),
        Some(w) => match w[1].as_str() {
            "fate" => paint::View::Fate,
            "oxygen" => paint::View::Oxygen,
            "clone" => paint::View::Clone,
            "stress" => paint::View::Stress,
            "vegf" => paint::View::Vegf,
            "flow" => paint::View::Flow,
            other => {
                return Err(anyhow!(
                    "--view {other}: expected fate, oxygen, clone, stress, vegf or flow"
                ));
            }
        },
    };
    let spec = config.load_network(config_path)?;
    let experiment = Experiment {
        label: config_path.unwrap_or("the default spheroid").to_owned(),
        scenario: config.scenario(config_path)?,
        spec,
    };
    let mut sides = vec![Side {
        experiment: experiment.clone(),
        name: "as it ran".to_owned(),
    }];
    if let Some(w) = args.windows(2).find(|w| w[0] == "--beside") {
        let path = w[1].as_str();
        let other =
            SimConfig::from_file(path).map_err(|e| anyhow!("loading config {path}: {e}"))?;
        let mut scenario = other.scenario(Some(path))?;
        scenario.seed = experiment.scenario.seed;
        if (scenario.width, scenario.height)
            != (experiment.scenario.width, experiment.scenario.height)
        {
            return Err(anyhow!(
                "--beside {path}: the two boxes must be the same size"
            ));
        }
        sides.push(Side {
            experiment: Experiment {
                label: path.to_owned(),
                scenario,
                spec: other.load_network(Some(path))?,
            },
            name: std::path::Path::new(path)
                .file_stem()
                .map_or_else(|| path.to_owned(), |s| s.to_string_lossy().into_owned()),
        });
    } else if let Some(w) = args.windows(2).find(|w| w[0] == "--undo") {
        let clones = w[1]
            .split(',')
            .map(|c| c.trim().parse::<u32>().map(CloneId))
            .collect::<Result<Vec<_>, _>>()
            .map_err(|e| anyhow!("--undo takes clone numbers separated by commas: {e}"))?;
        let at = match args.windows(2).find(|w| w[0] == "--undo-at") {
            Some(w) => w[1]
                .parse()
                .map_err(|e| anyhow!("--undo-at must be a step: {e}"))?,
            None => experiment
                .scenario
                .lesion
                .as_ref()
                .map(|l| l.at)
                .ok_or_else(|| anyhow!("--undo needs --undo-at <step> when there is no lesion"))?,
        };
        let names: Vec<String> = clones.iter().map(|c| c.0.to_string()).collect();
        let mut scenario = experiment.scenario.clone();
        scenario.revert = Some(Reversion { clones, at });
        sides.push(Side {
            experiment: Experiment {
                scenario,
                ..experiment.clone()
            },
            name: format!("clone {} undone at step {at}", names.join(", ")),
        });
    }

    App::new()
        .add_plugins(DefaultPlugins.set(WindowPlugin {
            primary_window: Some(Window {
                title: format!("cellweave · {}", experiment.label),
                ..default()
            }),
            ..default()
        }))
        .add_plugins(EguiPlugin::default())
        .insert_resource(ClearColor(paint::BACKGROUND))
        .insert_resource(runner::Pace::default())
        .insert_resource(runner::Divergence::default())
        .insert_resource(orbit::Orbit::framing(&experiment.scenario, sides.len()))
        .insert_resource(Sides(sides))
        .insert_resource(view)
        .insert_resource(experiment)
        .add_systems(Startup, setup)
        .add_systems(
            Update,
            (
                runner::keys,
                runner::collect,
                runner::compare,
                // Between the two: the simulation is back from the pool and
                // not yet sent off again, the one moment it can be read.
                paint::repaint,
                glow::retarget,
                relief::retarget,
                runner::dispatch,
                orbit::steer,
                orbit::spin,
                relief::toggle,
                relief::animate,
                glow::animate,
            )
                .chain(),
        )
        .add_systems(EguiPrimaryContextPass, panel::draw)
        .run();
    Ok(())
}

/// One run to show, and how the panel names it.
#[derive(Clone)]
pub struct Side {
    /// What it runs.
    pub experiment: Experiment,
    /// Its name in the panel.
    pub name: String,
}

/// The runs to show, left to right, until [`setup`] gives each its place.
#[derive(Resource)]
struct Sides(Vec<Side>);

/// Pixels of dark between two runs side by side.
const GAP: f32 = 24.0;

/// How bright the lamp is, in lux. Bright enough for the relief to read,
/// dim enough that the glow of a living cell still shows against it.
const LAMP: f32 = 450.0;

/// The camera, the lamp, and the surface the tissue is drawn on.
fn setup(
    mut commands: Commands,
    mut meshes: ResMut<Assets<Mesh>>,
    mut materials: ResMut<Assets<StandardMaterial>>,
    mut images: ResMut<Assets<Image>>,
    experiment: Res<Experiment>,
    sides: Res<Sides>,
    orbit: Res<orbit::Orbit>,
) {
    commands.spawn((
        Camera3d::default(),
        Camera {
            clear_color: ClearColorConfig::Custom(paint::BACKGROUND),
            ..default()
        },
        bevy::camera::Hdr,
        Tonemapping::AcesFitted,
        // Strong enough that a living cell glows into its membrane and the
        // vessel bleeds light into the tissue it feeds.
        Bloom {
            intensity: 0.3,
            ..Bloom::NATURAL
        },
        orbit.transform(),
    ));

    let (w, h) = (experiment.scenario.width, experiment.scenario.height);
    let size = w.max(h) as f32;

    // Low and to one side, so every dome throws a shadow across its
    // neighbour and the relief reads from any angle.
    commands.spawn((
        DirectionalLight {
            illuminance: LAMP,
            shadow_maps_enabled: true,
            color: Color::srgb(1.0, 0.95, 0.88),
            ..default()
        },
        CascadeShadowConfigBuilder {
            num_cascades: 2,
            maximum_distance: size * 4.0,
            first_cascade_far_bound: size * 1.6,
            ..default()
        }
        .build(),
        Transform::from_xyz(-size, size * 0.4, size * 0.25).looking_at(Vec3::ZERO, Vec3::Y),
    ));
    commands.insert_resource(GlobalAmbientLight {
        color: Color::srgb(0.45, 0.55, 1.0),
        brightness: 40.0,
        ..default()
    });

    // Side by side along X, centred on the origin as a whole.
    let n = sides.0.len() as f32;
    let span = w as f32 + GAP;
    for (index, side) in sides.0.iter().enumerate() {
        let x = (index as f32 - (n - 1.0) / 2.0) * span;
        let canvas = images.add(paint::canvas(w, h));
        let glow = images.add(glow::image(w, h));
        let material = materials.add(StandardMaterial {
            base_color_texture: Some(canvas.clone()),
            emissive: LinearRgba::rgb(paint::GLOW, paint::GLOW, paint::GLOW),
            emissive_texture: Some(glow.clone()),
            perceptual_roughness: 0.75,
            reflectance: 0.15,
            ..default()
        });
        let surface = meshes.add(relief::surface(w, h));
        commands.spawn((
            Mesh3d(surface.clone()),
            MeshMaterial3d(material.clone()),
            Transform::from_xyz(x, 0.0, 0.0),
            // The surface rises and falls every frame; the bounds it was
            // spawned with would cull it at the wrong moment.
            NoFrustumCulling,
            runner::Runner::new(side.experiment.clone(), side.name.clone(), index),
            paint::Canvas {
                image: canvas,
                material: material.clone(),
                seen: None,
            },
            glow::Glow::new(glow, material, w, h),
            relief::Relief::new(surface, w, h),
        ));
        // A slab under each tissue, so it reads as a dish with an edge
        // rather than a sheet floating in the dark.
        commands.spawn((
            Mesh3d(meshes.add(Cuboid::new(w as f32 + 6.0, 4.0, h as f32 + 6.0))),
            MeshMaterial3d(materials.add(StandardMaterial {
                base_color: Color::srgb(0.02, 0.025, 0.04),
                perceptual_roughness: 0.3,
                ..default()
            })),
            Transform::from_xyz(x, -2.05, 0.0),
        ));
    }
}
