//! A camera that circles the tissue: drag to turn it, scroll to come closer.

use bevy::input::mouse::{AccumulatedMouseMotion, AccumulatedMouseScroll, MouseScrollUnit};
use bevy::prelude::*;
use bevy_egui::EguiContexts;
use cellweave_engine::Scenario;

/// Radians turned per pixel dragged.
const TURN: f32 = 0.005;
/// Fraction of the distance covered per line scrolled.
const ZOOM: f32 = 0.1;
/// Radians per second the camera turns on its own, when asked to.
const SPIN: f32 = 0.12;
/// Never quite flat on the tissue, never past straight overhead.
const LOWEST: f32 = 0.15;
const HIGHEST: f32 = 1.55;

/// Where the camera sits, around the centre of the tissue.
#[derive(Resource, Debug, Clone, Copy)]
pub struct Orbit {
    /// Angle around the vertical axis.
    yaw: f32,
    /// Angle above the plane of the tissue.
    pitch: f32,
    /// Distance from the centre.
    distance: f32,
    /// Bounds on the distance, from the size of the box.
    nearest: f32,
    farthest: f32,
    /// Turning on its own.
    pub spinning: bool,
}

impl Orbit {
    /// Tilted enough to see the tissue as a surface, close enough that the
    /// whole box fills most of the window.
    pub fn framing(scenario: &Scenario) -> Self {
        let size = scenario.width.max(scenario.height) as f32;
        Self {
            yaw: 0.0,
            pitch: 0.8,
            distance: size * 1.3,
            nearest: size * 0.1,
            farthest: size * 4.0,
            spinning: false,
        }
    }

    /// The camera's transform for this position.
    pub fn transform(&self) -> Transform {
        let around = Vec3::new(
            self.yaw.sin() * self.pitch.cos(),
            self.pitch.sin(),
            self.yaw.cos() * self.pitch.cos(),
        );
        Transform::from_translation(around * self.distance).looking_at(Vec3::ZERO, Vec3::Y)
    }
}

/// Turn and zoom from the mouse, unless the pointer is over the panel.
pub fn steer(
    mut orbit: ResMut<Orbit>,
    buttons: Res<ButtonInput<MouseButton>>,
    motion: Res<AccumulatedMouseMotion>,
    scroll: Res<AccumulatedMouseScroll>,
    mut contexts: EguiContexts,
    mut camera: Query<&mut Transform, With<Camera3d>>,
) -> Result {
    let ctx = contexts.ctx_mut()?;
    if ctx.is_pointer_over_egui() || ctx.egui_wants_pointer_input() {
        return Ok(());
    }

    let mut moved = false;
    if buttons.pressed(MouseButton::Left) && motion.delta != Vec2::ZERO {
        // Taking hold of the camera stops it turning on its own.
        orbit.spinning = false;
        orbit.yaw -= motion.delta.x * TURN;
        orbit.pitch = (orbit.pitch + motion.delta.y * TURN).clamp(LOWEST, HIGHEST);
        moved = true;
    }
    let lines = match scroll.unit {
        MouseScrollUnit::Line => scroll.delta.y,
        MouseScrollUnit::Pixel => scroll.delta.y / 40.0,
    };
    if lines != 0.0 {
        orbit.distance =
            (orbit.distance * (1.0 - ZOOM * lines)).clamp(orbit.nearest, orbit.farthest);
        moved = true;
    }

    if moved {
        for mut transform in &mut camera {
            *transform = orbit.transform();
        }
    }
    Ok(())
}

/// A toggles a slow turn around the tissue; while on, the camera circles at
/// [`SPIN`] radians a second.
pub fn spin(
    time: Res<Time>,
    keys: Res<ButtonInput<KeyCode>>,
    mut orbit: ResMut<Orbit>,
    mut contexts: EguiContexts,
    mut camera: Query<&mut Transform, With<Camera3d>>,
) -> Result {
    if !contexts.ctx_mut()?.egui_wants_keyboard_input() && keys.just_pressed(KeyCode::KeyA) {
        orbit.spinning = !orbit.spinning;
    }
    if orbit.spinning {
        orbit.yaw += SPIN * time.delta_secs();
        for mut transform in &mut camera {
            *transform = orbit.transform();
        }
    }
    Ok(())
}
