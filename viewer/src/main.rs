//! cellweave-viewer: watch a simulation run.
//!
//! Usage:
//!
//! ```text
//! cellweave-viewer [--config <path>] [--view fate|oxygen|clone|stress]
//! ```
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
use cellweave_core::NetworkSpec;
use cellweave_engine::Scenario;
use cellweave_io::SimConfig;

/// What is being watched: the scenario and the hypothesis, kept so that a
/// restart starts the very same run again.
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

    let config = match config_path {
        Some(p) => SimConfig::from_file(p).map_err(|e| anyhow!("loading config {p}: {e}"))?,
        None => SimConfig::default(),
    };
    let view = match args.windows(2).find(|w| w[0] == "--view") {
        None => paint::View::default(),
        Some(w) => match w[1].as_str() {
            "fate" => paint::View::Fate,
            "oxygen" => paint::View::Oxygen,
            "clone" => paint::View::Clone,
            "stress" => paint::View::Stress,
            other => {
                return Err(anyhow!(
                    "--view {other}: expected fate, oxygen, clone or stress"
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
        .insert_resource(runner::Runner::new(&experiment))
        .insert_resource(orbit::Orbit::framing(&experiment.scenario))
        .insert_resource(view)
        .insert_resource(experiment)
        .add_systems(Startup, setup)
        .add_systems(
            Update,
            (
                runner::keys,
                runner::collect,
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
        Transform::IDENTITY,
        // The surface rises and falls every frame; the bounds it was spawned
        // with would cull it at the wrong moment.
        NoFrustumCulling,
    ));
    // A slab under the tissue, so it reads as a dish with an edge rather
    // than a sheet floating in the dark.
    commands.spawn((
        Mesh3d(meshes.add(Cuboid::new(w as f32 + 6.0, 4.0, h as f32 + 6.0))),
        MeshMaterial3d(materials.add(StandardMaterial {
            base_color: Color::srgb(0.02, 0.025, 0.04),
            perceptual_roughness: 0.3,
            ..default()
        })),
        Transform::from_xyz(0.0, -2.05, 0.0),
    ));

    commands.insert_resource(paint::Canvas {
        image: canvas,
        material: material.clone(),
    });
    commands.insert_resource(glow::Glow::new(glow, material, w, h));
    commands.insert_resource(relief::Relief::new(surface, w, h));
}
