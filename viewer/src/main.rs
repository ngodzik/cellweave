//! cellweave-viewer: watch a simulation run.
//!
//! Usage:
//!
//! ```text
//! cellweave-viewer [--config <path>]
//! ```
//!
//! It runs the same [`cellweave_engine::Simulation`] the command line runs,
//! under the same seed, so what it shows at a given step is what the command
//! line writes at that step. The tissue is drawn as a plane in HDR, each cell
//! lit from within by what it is doing, so that the structure a run produces,
//! a starving core or a cord around its vessel, is visible at a glance.
//!
//! Space pauses, R restarts, the mouse turns and zooms the camera.

#![deny(unsafe_code)]

mod orbit;
mod paint;
mod panel;
mod runner;

use anyhow::{Result, anyhow};
use bevy::core_pipeline::tonemapping::Tonemapping;
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
    let spec = config.load_network(config_path)?;
    let experiment = Experiment {
        label: config_path.unwrap_or("the default spheroid").to_owned(),
        scenario: config.scenario(),
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
        .insert_resource(paint::View::default())
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
                runner::dispatch,
                orbit::steer,
            )
                .chain(),
        )
        .add_systems(EguiPrimaryContextPass, panel::draw)
        .run();
    Ok(())
}

/// The camera, and the plane the tissue is painted on.
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
    let canvas = images.add(paint::canvas(w, h));
    // Lit only from within: black under any light, its colour all emission,
    // so that brightness is what the painter says it is and not where a lamp
    // happens to be.
    let material = materials.add(StandardMaterial {
        base_color: Color::BLACK,
        emissive: LinearRgba::rgb(paint::GLOW, paint::GLOW, paint::GLOW),
        emissive_texture: Some(canvas.clone()),
        perceptual_roughness: 1.0,
        reflectance: 0.0,
        ..default()
    });
    commands.spawn((
        Mesh3d(meshes.add(Plane3d::default().mesh().size(w as f32, h as f32))),
        MeshMaterial3d(material.clone()),
        Transform::IDENTITY,
    ));
    commands.insert_resource(paint::Canvas {
        image: canvas,
        material,
    });
}
