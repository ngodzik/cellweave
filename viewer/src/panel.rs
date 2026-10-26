//! The panel over the tissue: what the population is doing, and the controls.

use bevy::prelude::*;
use bevy_egui::{EguiContexts, egui};
use cellweave_core::OxygenSource;
use cellweave_engine::Measurements;
use egui_plot::{Legend, Line, Plot, PlotPoints};

use crate::Experiment;
use crate::paint::{self, View};
use crate::runner::Runner;

/// Draw the panel, and apply what was clicked in it.
pub fn draw(
    mut contexts: EguiContexts,
    mut runner: ResMut<Runner>,
    mut view: ResMut<View>,
    experiment: Res<Experiment>,
) -> Result {
    let ctx = contexts.ctx_mut()?;
    egui::Window::new("cellweave")
        .default_pos([12.0, 12.0])
        .default_width(300.0)
        .resizable(false)
        .show(ctx, |ui| {
            describe(ui, &experiment);
            ui.separator();
            census(ui, &runner);
            ui.separator();
            controls(ui, &mut runner, &experiment);
            ui.separator();
            colouring(ui, &mut view, &experiment);
            ui.separator();
            plot(ui, &runner);
            if let Some(why) = &runner.stopped {
                ui.separator();
                ui.colored_label(colour(paint::QUIESCENT), why);
            }
            ui.separator();
            ui.weak("Space pause · N step · R restart · drag to turn · scroll to zoom");
        });
    Ok(())
}

fn describe(ui: &mut egui::Ui, experiment: &Experiment) {
    let s = &experiment.scenario;
    let source = match s.oxygen_source {
        OxygenSource::Bath => "a bath around the tissue".to_owned(),
        OxygenSource::Vessel { radius } => format!("a vessel of radius {radius}"),
    };
    ui.strong(&experiment.label);
    ui.label(format!(
        "{}×{} pixels, oxygen from {source}, seed {}",
        s.width, s.height, s.seed
    ));
}

fn census(ui: &mut egui::Ui, runner: &Runner) {
    let Some(c) = runner.latest() else {
        ui.label(if runner.settling() {
            "Solving the oxygen field for the first time…"
        } else {
            "No step taken yet."
        });
        return;
    };
    egui::Grid::new("census").num_columns(2).show(ui, |ui| {
        let mut row = |name: &str, value: String| {
            ui.label(name);
            ui.monospace(value);
            ui.end_row();
        };
        row("MCS", c.mcs.to_string());
        row("living", format!("{} ({} cycling)", c.living(), c.cycling));
        row("necrotic", c.necrotic.to_string());
        row("born · gone", format!("{} · {}", c.born, c.gone));
        row("viable depth", format!("{:.1} px", c.viable_depth));
        row("necrotic radius", format!("{:.1} px", c.necrotic_radius));
        row("speed", format!("{:.0} steps/s", runner.rate));
    });
}

fn controls(ui: &mut egui::Ui, runner: &mut Runner, experiment: &Experiment) {
    ui.horizontal(|ui| {
        let label = if runner.paused { "Run" } else { "Pause" };
        if ui.button(label).clicked() {
            runner.paused = !runner.paused;
        }
        if ui.button("Step").clicked() {
            runner.paused = true;
            runner.step_once = true;
        }
        if ui.button("Restart").clicked() {
            runner.restart(experiment);
        }
    });
    ui.add(egui::Slider::new(&mut runner.steps_per_batch, 1..=25).text("steps per batch"));
}

fn colouring(ui: &mut egui::Ui, view: &mut View, experiment: &Experiment) {
    let mut chosen = *view;
    ui.horizontal(|ui| {
        ui.label("Colour by");
        ui.selectable_value(&mut chosen, View::Fate, "fate");
        ui.selectable_value(&mut chosen, View::Oxygen, "oxygen");
        let protein = match chosen {
            View::Protein(i) => i,
            _ => 0,
        };
        ui.selectable_value(&mut chosen, View::Protein(protein), "protein");
    });

    match chosen {
        View::Fate => {
            ui.horizontal_wrapped(|ui| {
                swatch(ui, paint::CYCLING, "cycling");
                swatch(ui, paint::QUIESCENT, "quiescent");
                swatch(ui, paint::NECROTIC, "necrotic");
                swatch(ui, paint::VESSEL, "vessel");
            });
        }
        View::Protein(i) => {
            let names: Vec<&str> = experiment
                .spec
                .nodes
                .iter()
                .map(|n| n.name.as_str())
                .collect();
            let mut index = i.min(names.len().saturating_sub(1));
            egui::ComboBox::from_label("protein")
                .selected_text(names.get(index).copied().unwrap_or("none"))
                .show_ui(ui, |ui| {
                    for (k, name) in names.iter().enumerate() {
                        ui.selectable_value(&mut index, k, *name);
                    }
                });
            chosen = View::Protein(index);
            gradient(ui);
        }
        View::Oxygen => gradient(ui),
    }

    // Assigned only on a real change, so the canvas is not repainted every
    // frame the panel is drawn.
    if chosen != *view {
        *view = chosen;
    }
}

/// The population, then the shape: how many cells, and how deep the living
/// reach.
fn plot(ui: &mut egui::Ui, runner: &Runner) {
    let series = |pick: fn(&Measurements) -> f64| -> PlotPoints {
        runner
            .history
            .iter()
            .map(|m| [m.mcs as f64, pick(m)])
            .collect()
    };
    let frame = |id: &str| {
        Plot::new(id.to_owned())
            .height(100.0)
            .legend(Legend::default())
            .allow_drag(false)
            .allow_zoom(false)
            .allow_scroll(false)
            .allow_boxed_zoom(false)
            .include_y(0.0)
    };
    frame("population").show(ui, |plot| {
        plot.line(Line::new("living", series(|m| m.living() as f64)).color(colour(paint::CYCLING)));
        plot.line(
            Line::new("necrotic", series(|m| m.necrotic as f64))
                .color(egui::Color32::from_gray(150)),
        );
    });
    frame("shape").show(ui, |plot| {
        plot.line(
            Line::new("viable depth", series(|m| m.viable_depth)).color(colour(paint::QUIESCENT)),
        );
        plot.line(
            Line::new("necrotic radius", series(|m| m.necrotic_radius))
                .color(egui::Color32::from_gray(150)),
        );
    });
}

/// A small square of colour with its meaning.
fn swatch(ui: &mut egui::Ui, fill: LinearRgba, name: &str) {
    let (rect, _) = ui.allocate_exact_size(egui::vec2(12.0, 12.0), egui::Sense::hover());
    ui.painter().rect_filled(rect, 2.0, colour(fill));
    ui.label(name);
}

/// The ramp protein levels and oxygen are drawn with, from zero to one.
fn gradient(ui: &mut egui::Ui) {
    ui.horizontal(|ui| {
        ui.label("0");
        let (rect, _) = ui.allocate_exact_size(egui::vec2(160.0, 10.0), egui::Sense::hover());
        let steps = 32;
        for k in 0..steps {
            let t = k as f32 / (steps - 1) as f32;
            let x0 = rect.left() + rect.width() * k as f32 / steps as f32;
            let x1 = rect.left() + rect.width() * (k + 1) as f32 / steps as f32;
            let slice = egui::Rect::from_x_y_ranges(x0..=x1, rect.y_range());
            ui.painter().rect_filled(slice, 0.0, colour(paint::ramp(t)));
        }
        ui.label("1");
    });
}

fn colour(c: LinearRgba) -> egui::Color32 {
    let [r, g, b, _] = paint::bytes(c);
    egui::Color32::from_rgb(r, g, b)
}
