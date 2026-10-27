//! Steps the simulation in the task pool, so the window keeps turning
//! smoothly while the field settles, which is most of what a step costs.
//!
//! The simulation itself moves into the task and comes back out with it: it is
//! owned by exactly one place at a time, and nothing is shared or locked.

use std::time::Instant;

use bevy::prelude::*;
use bevy::tasks::futures::check_ready;
use bevy::tasks::{AsyncComputeTaskPool, Task};
use bevy_egui::EguiContexts;
use cellweave_engine::Simulation;
use cellweave_engine::simulation::{Configured, Running, StepReport};

use crate::Experiment;

/// Where the simulation is.
enum Stage {
    /// Built, and waiting to have its field solved in the pool.
    Configured(Simulation<Configured>),
    /// Back from the pool, ready for more steps, and what is painted.
    Idle(Simulation<Running>),
    /// In the pool.
    Busy(Task<Outcome>),
    /// Lost to an error; only a restart brings it back.
    Failed,
}

/// What the pool is asked to do.
enum Work {
    Start(Simulation<Configured>),
    Steps(Simulation<Running>, u32),
}

/// What comes back from the pool.
struct Outcome {
    sim: Result<Simulation<Running>, String>,
    reports: Vec<StepReport>,
    began: Instant,
}

/// Run one batch of work. Stops at the first step that reaches the edge of
/// the box, since every step after that would be an artefact of its size.
fn work(job: Work) -> Outcome {
    let began = Instant::now();
    let mut reports = Vec::new();
    let sim = match job {
        Work::Start(sim) => sim.start().map(|(sim, _)| sim),
        Work::Steps(mut sim, steps) => (|| {
            for _ in 0..steps {
                let report = sim.step()?;
                reports.push(report);
                if report.reached_edge {
                    break;
                }
            }
            Ok(sim)
        })(),
    };
    Outcome {
        sim: sim.map_err(|e: cellweave_core::CellweaveError| e.to_string()),
        reports,
        began,
    }
}

/// The population at one step, for the panel and its plot.
#[derive(Debug, Clone, Copy)]
pub struct Census {
    /// The step it was taken at.
    pub mcs: u64,
    /// Alive, cycling or quiescent.
    pub living: usize,
    /// Dead and not yet resorbed.
    pub necrotic: usize,
    /// Ever born.
    pub born: u64,
    /// Resorbed to nothing.
    pub gone: u64,
}

/// The simulation, wherever it is, and what the panel asks of it.
#[derive(Resource)]
pub struct Runner {
    stage: Stage,
    /// The canvas is out of date.
    fresh: bool,
    /// Nothing is stepped while set, except one step when asked.
    pub paused: bool,
    /// Take exactly one step, then wait.
    pub step_once: bool,
    /// Steps per batch sent to the pool. More is faster to watch and coarser
    /// to see.
    pub steps_per_batch: u32,
    /// Why the run stopped, if it did.
    pub stopped: Option<String>,
    /// One census per batch, oldest first.
    pub history: Vec<Census>,
    /// Steps per second, smoothed.
    pub rate: f64,
}

impl Runner {
    /// A runner holding a fresh simulation of this experiment.
    pub fn new(experiment: &Experiment) -> Self {
        let mut runner = Self {
            stage: Stage::Failed,
            fresh: true,
            paused: false,
            step_once: false,
            steps_per_batch: 2,
            stopped: None,
            history: Vec::new(),
            rate: 0.0,
        };
        runner.restart(experiment);
        runner
    }

    /// Throw the current run away and start the same one again. A batch in
    /// the pool is dropped, which cancels it.
    pub fn restart(&mut self, experiment: &Experiment) {
        self.history.clear();
        self.rate = 0.0;
        self.fresh = true;
        match Simulation::new(experiment.scenario.clone(), &experiment.spec) {
            Ok(sim) => {
                self.stage = Stage::Configured(sim);
                self.stopped = None;
            }
            Err(e) => {
                self.stage = Stage::Failed;
                self.stopped = Some(format!("The run could not be built: {e}"));
            }
        }
    }

    /// The simulation, when it is here and not in the pool.
    pub fn idle(&self) -> Option<&Simulation<Running>> {
        match &self.stage {
            Stage::Idle(sim) => Some(sim),
            _ => None,
        }
    }

    /// Whether the field is still being solved for the first time.
    pub fn settling(&self) -> bool {
        matches!(self.stage, Stage::Configured(_))
            || (matches!(self.stage, Stage::Busy(_)) && self.history.is_empty())
    }

    /// Whether the canvas needs painting, clearing the flag.
    pub fn take_fresh(&mut self) -> bool {
        std::mem::take(&mut self.fresh)
    }

    /// The latest census, if a step has been taken.
    pub fn latest(&self) -> Option<&Census> {
        self.history.last()
    }
}

/// Space pauses, N steps once, R restarts, unless egui is taking the keys.
pub fn keys(
    mut runner: ResMut<Runner>,
    keys: Res<ButtonInput<KeyCode>>,
    experiment: Res<Experiment>,
    mut contexts: EguiContexts,
) -> Result {
    if contexts.ctx_mut()?.egui_wants_keyboard_input() {
        return Ok(());
    }
    if keys.just_pressed(KeyCode::Space) {
        runner.paused = !runner.paused;
    }
    if keys.just_pressed(KeyCode::KeyN) {
        runner.paused = true;
        runner.step_once = true;
    }
    if keys.just_pressed(KeyCode::KeyR) {
        runner.restart(&experiment);
    }
    Ok(())
}

/// Take back a batch the pool has finished, if there is one.
pub fn collect(mut runner: ResMut<Runner>) {
    let runner = &mut *runner;
    let Stage::Busy(task) = &mut runner.stage else {
        return;
    };
    let Some(outcome) = check_ready(task) else {
        return;
    };

    if let Some(edge) = outcome.reports.iter().find(|r| r.reached_edge) {
        runner.stopped = Some(format!(
            "The tissue reached the edge of the box at MCS {}. Everything after would be an \
             artefact of the box's size; a larger grid goes further.",
            edge.mcs
        ));
    }
    if !outcome.reports.is_empty() {
        let seconds = outcome.began.elapsed().as_secs_f64().max(1e-6);
        let now = outcome.reports.len() as f64 / seconds;
        runner.rate = if runner.rate == 0.0 {
            now
        } else {
            0.8 * runner.rate + 0.2 * now
        };
    }

    match outcome.sim {
        Ok(sim) => {
            if let Some(last) = outcome.reports.last() {
                runner.history.push(Census {
                    mcs: last.mcs,
                    living: sim.living(),
                    necrotic: sim.necrotic(),
                    born: sim.born(),
                    gone: sim.gone(),
                });
            }
            runner.stage = Stage::Idle(sim);
        }
        Err(e) => {
            runner.stage = Stage::Failed;
            runner.stopped = Some(format!("The run failed: {e}"));
        }
    }
    runner.fresh = true;
}

/// Hand the simulation to the pool, if there is something for it to do.
pub fn dispatch(mut runner: ResMut<Runner>) {
    let runner = &mut *runner;
    let job = match std::mem::replace(&mut runner.stage, Stage::Failed) {
        Stage::Configured(sim) => Work::Start(sim),
        Stage::Idle(sim) if runner.stopped.is_none() && (!runner.paused || runner.step_once) => {
            let steps = if runner.step_once {
                1
            } else {
                runner.steps_per_batch
            };
            runner.step_once = false;
            Work::Steps(sim, steps)
        }
        other => {
            runner.stage = other;
            return;
        }
    };
    let task = AsyncComputeTaskPool::get().spawn(async move { work(job) });
    runner.stage = Stage::Busy(task);
}
