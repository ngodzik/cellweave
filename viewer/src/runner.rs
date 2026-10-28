//! Steps the simulation in the task pool, so the window keeps turning
//! smoothly while the field settles, which is most of what a step costs.
//!
//! The simulation itself moves into the task and comes back out with it: it is
//! owned by exactly one place at a time, and nothing is shared or locked.
//!
//! Two runs can be shown side by side, each with a runner of its own. They are
//! stepped together, batch for batch, so that both always show the same step.

use std::time::Instant;

use bevy::prelude::*;
use bevy::tasks::futures::check_ready;
use bevy::tasks::{AsyncComputeTaskPool, Task};
use bevy_egui::EguiContexts;
use cellweave_engine::simulation::{Configured, Running, StepReport};
use cellweave_engine::{CloneRecord, Measurements, Simulation};

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

/// Entries of history kept at most. Past this, every other one is dropped,
/// so a long run keeps its whole span at a coarser grain, and what the panel
/// draws each frame stays the same size however long the window is open.
const HISTORY_LIMIT: usize = 2000;

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

/// How the runs are stepped: shared by every run on screen, so that they
/// pause, step and go together.
#[derive(Resource, Debug)]
pub struct Pace {
    /// Nothing is stepped while set, except one step when asked.
    pub paused: bool,
    /// Take exactly one step, then wait.
    pub step_once: bool,
    /// Steps per batch sent to the pool. More is faster to watch and coarser
    /// to see.
    pub steps_per_batch: u32,
}

impl Default for Pace {
    fn default() -> Self {
        Self {
            paused: false,
            step_once: false,
            steps_per_batch: 2,
        }
    }
}

/// Where two runs shown side by side stopped being the same run.
#[derive(Resource, Debug, Default)]
pub struct Divergence {
    /// The step of the first batch whose measurements differed, once they
    /// have; `None` while they are identical, or with one run.
    pub parted_at: Option<u64>,
}

/// One run on screen: the simulation, wherever it is, and what the panel
/// asks of it.
#[derive(Component)]
pub struct Runner {
    /// What it runs, kept so that a restart starts the very same run again.
    pub experiment: Experiment,
    /// How the panel names it.
    pub name: String,
    /// Its place on screen, from the left; the first is the one the panel
    /// describes in full.
    pub index: usize,
    stage: Stage,
    /// Counts every change to what the simulation shows: each batch back from
    /// the pool, each restart. Anything drawn from the simulation remembers
    /// the generation it drew, and draws again when this moves on.
    generation: u64,
    /// The generation the last restart began.
    restarts_at: u64,
    /// Why the run stopped, if it did.
    pub stopped: Option<String>,
    /// The tissue measured once per batch, oldest first.
    pub history: Vec<Measurements>,
    /// Living cells per clone, once per batch, alongside `history`.
    pub clone_history: Vec<Vec<usize>>,
    /// Every clone founded so far, as of the last batch, so the panel can
    /// name them while the simulation is away in the pool.
    pub clones: Vec<CloneRecord>,
    /// Steps per second, smoothed.
    pub rate: f64,
}

impl Runner {
    /// A runner holding a fresh simulation of this experiment.
    pub fn new(experiment: Experiment, name: String, index: usize) -> Self {
        let mut runner = Self {
            experiment,
            name,
            index,
            stage: Stage::Failed,
            generation: 0,
            restarts_at: 0,
            stopped: None,
            history: Vec::new(),
            clone_history: Vec::new(),
            clones: Vec::new(),
            rate: 0.0,
        };
        runner.restart();
        runner
    }

    /// Throw the current run away and start the same one again. A batch in
    /// the pool is dropped, which cancels it.
    pub fn restart(&mut self) {
        self.history.clear();
        self.clone_history.clear();
        self.clones.clear();
        self.rate = 0.0;
        self.generation += 1;
        self.restarts_at = self.generation;
        match Simulation::new(self.experiment.scenario.clone(), &self.experiment.spec) {
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

    /// Whether a batch is away in the pool.
    pub fn busy(&self) -> bool {
        matches!(self.stage, Stage::Busy(_))
    }

    /// Whether the field is still being solved for the first time.
    pub fn settling(&self) -> bool {
        matches!(self.stage, Stage::Configured(_))
            || (matches!(self.stage, Stage::Busy(_)) && self.history.is_empty())
    }

    /// Which state of the simulation is the current one; see the field.
    pub fn generation(&self) -> u64 {
        self.generation
    }

    /// Whether the run was restarted since `seen`, rather than stepped: a
    /// census of the population then starts from nothing.
    pub fn restarted_since(&self, seen: u64) -> bool {
        self.restarts_at > seen
    }

    /// Halve the history, in both series together so they stay aligned.
    fn thin_history(&mut self) {
        thin(&mut self.history);
        thin(&mut self.clone_history);
    }

    /// The latest measurements, if a step has been taken.
    pub fn latest(&self) -> Option<&Measurements> {
        self.history.last()
    }
}

/// Keep every other entry of a series, and always the latest.
fn thin<T>(series: &mut Vec<T>) {
    let last = series.len().saturating_sub(1);
    let mut i = 0;
    series.retain(|_| {
        let keep = i % 2 == 0 || i == last;
        i += 1;
        keep
    });
}

/// Space pauses, N steps once, R restarts every run, unless egui is taking
/// the keys.
pub fn keys(
    mut runners: Query<&mut Runner>,
    mut pace: ResMut<Pace>,
    keys: Res<ButtonInput<KeyCode>>,
    mut contexts: EguiContexts,
) -> Result {
    if contexts.ctx_mut()?.egui_wants_keyboard_input() {
        return Ok(());
    }
    if keys.just_pressed(KeyCode::Space) {
        pace.paused = !pace.paused;
    }
    if keys.just_pressed(KeyCode::KeyN) {
        pace.paused = true;
        pace.step_once = true;
    }
    if keys.just_pressed(KeyCode::KeyR) {
        for mut runner in &mut runners {
            runner.restart();
        }
    }
    Ok(())
}

/// Take back every batch the pool has finished.
pub fn collect(mut runners: Query<&mut Runner>) {
    for mut runner in &mut runners {
        collect_one(&mut runner);
    }
}

/// Note where two runs side by side first differ: the step of the first pair
/// of batches whose measurements are not the same.
pub fn compare(runners: Query<&Runner>, mut divergence: ResMut<Divergence>) {
    let mut sides: Vec<&Runner> = runners.iter().collect();
    sides.sort_by_key(|r| r.index);
    let [a, b] = sides.as_slice() else {
        divergence.parted_at = None;
        return;
    };
    if a.history.is_empty() || b.history.is_empty() {
        divergence.parted_at = None;
        return;
    }
    if divergence.parted_at.is_none()
        && let (Some(x), Some(y)) = (a.latest(), b.latest())
        && x.mcs == y.mcs
        && x != y
    {
        divergence.parted_at = Some(x.mcs);
    }
}

fn collect_one(runner: &mut Runner) {
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
            if !outcome.reports.is_empty() {
                runner.history.push(sim.measure());
                runner.clone_history.push(sim.clone_sizes());
                if runner.history.len() > HISTORY_LIMIT {
                    runner.thin_history();
                }
                if runner.clones.len() != sim.clones().len() {
                    runner.clones = sim.clones().to_vec();
                }
            }
            runner.stage = Stage::Idle(sim);
        }
        Err(e) => {
            runner.stage = Stage::Failed;
            runner.stopped = Some(format!("The run failed: {e}"));
        }
    }
    runner.generation += 1;
}

/// Hand every run to the pool, if there is something for it to do, and only
/// once none is still away: runs side by side move batch for batch, and so
/// always show the same step.
pub fn dispatch(mut runners: Query<&mut Runner>, mut pace: ResMut<Pace>) {
    if runners.iter().any(Runner::busy) {
        return;
    }
    let go = !pace.paused || pace.step_once;
    let steps = if pace.step_once {
        1
    } else {
        pace.steps_per_batch
    };
    let mut sent = false;
    for mut runner in &mut runners {
        let runner = &mut *runner;
        let job = match std::mem::replace(&mut runner.stage, Stage::Failed) {
            Stage::Configured(sim) => Work::Start(sim),
            Stage::Idle(sim) if runner.stopped.is_none() && go => {
                sent = true;
                Work::Steps(sim, steps)
            }
            other => {
                runner.stage = other;
                continue;
            }
        };
        let task = AsyncComputeTaskPool::get().spawn(async move { work(job) });
        runner.stage = Stage::Busy(task);
    }
    if sent {
        pace.step_once = false;
    }
}
