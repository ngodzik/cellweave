//! Saving what the window shows as the run goes, to make an animation of it.
//!
//! The frames are taken by the renderer from the image it draws, not grabbed
//! from the screen, so whatever lies over the window on screen never gets
//! into them.

use std::path::PathBuf;

use bevy::prelude::*;
use bevy::render::view::screenshot::{Screenshot, save_to_disk};

use crate::runner::Runner;

/// Frames the renderer is given to save the last one before the viewer
/// quits: a screenshot is captured on a later frame than the one it is asked
/// for.
const FRAMES_TO_FINISH: u8 = 10;

/// Where frames go, and when to stop.
#[derive(Resource, Debug)]
pub struct Recording {
    /// The directory every frame is written to, as `frame_<step>.png`.
    dir: PathBuf,
    /// The step after which the viewer quits, once its frame is taken.
    until: Option<u64>,
    /// The last step a frame was taken of.
    last: Option<u64>,
    /// Frames left before quitting, once the last one has been asked for.
    finishing: Option<u8>,
}

impl Recording {
    /// Frames into `dir`, quitting after step `until` if there is one.
    pub fn new(dir: PathBuf, until: Option<u64>) -> Self {
        Self {
            dir,
            until,
            last: None,
            finishing: None,
        }
    }
}

/// Take a frame each time the first run shows a step it had not shown yet,
/// and quit a few frames after the last one asked for, or after the run
/// stops.
pub fn capture(
    mut commands: Commands,
    recording: Option<ResMut<Recording>>,
    runners: Query<&Runner>,
    mut exit: MessageWriter<AppExit>,
) {
    let Some(mut recording) = recording else {
        return;
    };
    if let Some(left) = recording.finishing {
        if left == 0 {
            exit.write(AppExit::Success);
        } else {
            recording.finishing = Some(left - 1);
        }
        return;
    }
    let Some(first) = runners.iter().find(|r| r.index == 0) else {
        return;
    };
    let Some(step) = first.history.last().map(|m| m.mcs) else {
        return;
    };
    if recording.last.is_some_and(|last| step <= last) {
        return;
    }
    recording.last = Some(step);
    let path = recording.dir.join(format!("frame_{step:06}.png"));
    commands
        .spawn(Screenshot::primary_window())
        .observe(save_to_disk(path));
    // A run that stops, at the edge of its box or failing, takes no more
    // steps: its last frame is the last there will be.
    if recording.until.is_some_and(|until| step >= until) || first.stopped.is_some() {
        recording.finishing = Some(FRAMES_TO_FINISH);
    }
}
