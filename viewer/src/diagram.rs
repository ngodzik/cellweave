//! A cell's network drawn as a diagram: what it reads on the left, each
//! protein in a column by how many relays it sits from what it reads, and
//! every term that reads something as an arrow, as thick as it pushes.
//!
//! The layout is worked out from the network itself, not drawn by hand, so
//! that any network a configuration names, and any cell whose mutations
//! rewired it, is drawn as it is.

use bevy_egui::egui::{self, Color32, Pos2, Rect, Sense, Shape, Stroke, Vec2};
use cellweave_core::network::INPUT_NAMES;

use crate::inspect::Link;
use crate::paint;

/// Room per column and per row, and the size of a protein's box.
const COLUMN: f32 = 90.0;
const ROW: f32 = 60.0;
const BOX: Vec2 = Vec2::new(76.0, 36.0);
/// Room around the whole, for the arrows that loop back under a row.
const MARGIN: f32 = 18.0;

/// A rate of change, per unit of network time, at which an arrow is drawn
/// at its thickest.
const FULL_PUSH: f64 = 0.15;

/// Activation, inhibition and damage: green, red, and a darker red.
const ACTIVATES: Color32 = Color32::from_rgb(90, 205, 150);
const INHIBITS: Color32 = Color32::from_rgb(240, 95, 90);
const DAMAGES: Color32 = Color32::from_rgb(190, 60, 120);

/// Where everything goes, in the diagram's own coordinates.
pub struct Layout {
    /// The inputs the network reads, by name, in the engine's order.
    pub inputs: Vec<&'static str>,
    /// The centre of each input, then of each node, in that order.
    pub centres: Vec<Pos2>,
    /// The column of each, the same order.
    columns: Vec<usize>,
    /// The size of the whole.
    pub size: Vec2,
}

impl Layout {
    /// Lay out `nodes` proteins and what `links` say they read.
    pub fn of(nodes: usize, names: &[String], links: &[Link]) -> Self {
        let inputs: Vec<&'static str> = INPUT_NAMES
            .iter()
            .copied()
            .filter(|input| links.iter().any(|l| l.source.as_deref() == Some(*input)))
            .collect();
        let vertices = inputs.len() + nodes;
        let index = |name: &str| {
            inputs.iter().position(|i| *i == name).or_else(|| {
                names
                    .iter()
                    .position(|n| n == name)
                    .map(|n| inputs.len() + n)
            })
        };
        let mut edges: Vec<(usize, usize)> = links
            .iter()
            .filter_map(|l| Some((index(l.source.as_deref()?)?, inputs.len() + l.node)))
            .filter(|(from, to)| from != to)
            .collect();
        edges.sort_unstable();
        edges.dedup();

        let columns = columns(vertices, inputs.len(), &edges);
        let rows = rows(&columns, &edges);
        let tallest = rows.iter().copied().fold(0, usize::max) + 1;
        let mut per_column = vec![0usize; columns.iter().copied().max().unwrap_or(0) + 1];
        for &c in &columns {
            per_column[c] += 1;
        }
        let centres = (0..vertices)
            .map(|v| {
                // Each column centred on the tallest one.
                let offset = (tallest - per_column[columns[v]]) as f32 / 2.0;
                Pos2::new(
                    MARGIN + BOX.x / 2.0 + columns[v] as f32 * COLUMN,
                    MARGIN + BOX.y / 2.0 + (rows[v] as f32 + offset) * ROW,
                )
            })
            .collect();
        Self {
            inputs,
            centres,
            columns,
            size: Vec2::new(
                2.0 * MARGIN + BOX.x + (per_column.len() - 1) as f32 * COLUMN,
                2.0 * MARGIN + BOX.y + (tallest - 1) as f32 * ROW + 16.0,
            ),
        }
    }

    /// The vertex a term's source is, if it reads one.
    fn source_of(&self, link: &Link, names: &[String]) -> Option<usize> {
        let name = link.source.as_deref()?;
        self.inputs.iter().position(|i| *i == name).or_else(|| {
            names
                .iter()
                .position(|n| n == name)
                .map(|n| self.inputs.len() + n)
        })
    }
}

/// The column of every vertex: inputs in the first, and each node one past
/// the furthest of what it reads, once the links that close a loop back are
/// set aside. A node that reads nothing goes in the second.
fn columns(vertices: usize, inputs: usize, edges: &[(usize, usize)]) -> Vec<usize> {
    let mut after: Vec<Vec<usize>> = vec![Vec::new(); vertices];
    for &(from, to) in edges {
        after[from].push(to);
    }
    // A depth first walk from the inputs, then from whatever it missed: an
    // edge to a vertex still on the walk's path is one that loops back.
    let mut state = vec![0u8; vertices];
    let mut back = Vec::new();
    fn walk(v: usize, after: &[Vec<usize>], state: &mut [u8], back: &mut Vec<(usize, usize)>) {
        state[v] = 1;
        for &next in &after[v] {
            match state[next] {
                0 => walk(next, after, state, back),
                1 => back.push((v, next)),
                _ => {}
            }
        }
        state[v] = 2;
    }
    for v in 0..vertices {
        if state[v] == 0 {
            walk(v, &after, &mut state, &mut back);
        }
    }
    let forward: Vec<(usize, usize)> = edges
        .iter()
        .copied()
        .filter(|e| !back.contains(e))
        .collect();
    let mut column: Vec<usize> = (0..vertices).map(|v| usize::from(v >= inputs)).collect();
    // The longest path, by relaxing every edge as many times as there are
    // vertices: a graph this small needs no cleverer order.
    for _ in 0..vertices {
        for &(from, to) in &forward {
            if to >= inputs {
                column[to] = column[to].max(column[from] + 1);
            }
        }
    }
    column
}

/// The row of every vertex within its column, ordered so that each sits
/// near what it is linked to: a few sweeps, each sorting a column by the
/// mean row of its neighbours.
fn rows(columns: &[usize], edges: &[(usize, usize)]) -> Vec<usize> {
    let n = columns.len();
    let mut row = vec![0.0f32; n];
    let count = columns.iter().copied().max().unwrap_or(0) + 1;
    let mut members: Vec<Vec<usize>> = vec![Vec::new(); count];
    for v in 0..n {
        members[columns[v]].push(v);
    }
    for list in &members {
        for (r, &v) in list.iter().enumerate() {
            row[v] = r as f32;
        }
    }
    for _ in 0..6 {
        for list in &mut members {
            let mean = |v: usize| {
                let (sum, k) = edges
                    .iter()
                    .filter_map(|&(a, b)| match (a == v, b == v) {
                        (true, false) => Some(row[b]),
                        (false, true) => Some(row[a]),
                        _ => None,
                    })
                    .fold((0.0, 0), |(s, k), r| (s + r, k + 1));
                if k == 0 { row[v] } else { sum / k as f32 }
            };
            let mut keyed: Vec<(f32, usize)> = list.iter().map(|&v| (mean(v), v)).collect();
            keyed.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
            *list = keyed.into_iter().map(|(_, v)| v).collect();
            for (r, &v) in list.iter().enumerate() {
                row[v] = r as f32;
            }
        }
    }
    row.into_iter().map(|r| r as usize).collect()
}

/// What a click on the diagram asked for.
pub enum Clicked {
    /// Nothing.
    Nothing,
    /// The protein with this index.
    Node(usize),
}

/// Draw the diagram, scaled to fit `width`, and say what was clicked.
///
/// `levels` are the proteins' levels, `reads` the inputs' values by name,
/// `badges` a short role label per protein, and `chosen` the protein the
/// tissue is coloured by, if one is.
#[allow(clippy::too_many_arguments)]
pub fn draw(
    ui: &mut egui::Ui,
    layout: &Layout,
    names: &[String],
    levels: &[f64],
    links: &[Link],
    reads: &dyn Fn(&str) -> f64,
    badges: &[Option<(String, bool)>],
    chosen: Option<usize>,
    explain: &dyn Fn(usize) -> String,
) -> Clicked {
    let scale = (ui.available_width() / layout.size.x).min(1.0);
    let (whole, _) = ui.allocate_exact_size(layout.size * scale, Sense::hover());
    let painter = ui.painter_at(whole);
    let at = |p: Pos2| whole.min + (p.to_vec2() * scale);
    let half = BOX * scale / 2.0;
    let boxed = |v: usize| Rect::from_center_size(at(layout.centres[v]), BOX * scale);
    let inputs = layout.inputs.len();

    // Arrows first, under the boxes.
    for link in links {
        let Some(from) = layout.source_of(link, names) else {
            continue;
        };
        let to = inputs + link.node;
        if from == to {
            continue;
        }
        let strength = ((link.now.abs() / FULL_PUSH) as f32).min(1.0);
        let (colour, dashed, bar) = match link.kind {
            "activated_by" | "produced_by" => (ACTIVATES, false, false),
            "activated_by_absence_of" => (ACTIVATES, true, false),
            "inhibited_by" => (INHIBITS, false, true),
            _ => (DAMAGES, true, true),
        };
        let colour = colour.gamma_multiply(0.35 + 0.65 * strength);
        let stroke = Stroke::new((1.0 + 3.5 * strength) * scale.max(0.6), colour);
        let (a, b) = (at(layout.centres[from]), at(layout.centres[to]));
        let forward = layout.columns[to] > layout.columns[from];
        let points = if forward {
            let start = a + Vec2::new(half.x, 0.0);
            let end = b - Vec2::new(half.x + 4.0 * scale, 0.0);
            let bend = (end.x - start.x) / 2.0;
            [
                start,
                start + Vec2::new(bend, 0.0),
                end - Vec2::new(bend, 0.0),
                end,
            ]
        } else {
            // A loop back, under both boxes.
            let start = a + Vec2::new(0.0, half.y);
            let end = b + Vec2::new(0.0, half.y + 4.0 * scale);
            let drop = Vec2::new(0.0, 0.9 * ROW * scale);
            [start, start + drop, end + drop, end]
        };
        let curve = egui::epaint::CubicBezierShape::from_points_stroke(
            points,
            false,
            Color32::TRANSPARENT,
            stroke,
        )
        .flatten(Some(0.5));
        if dashed {
            painter.extend(Shape::dashed_line(&curve, stroke, 5.0 * scale, 4.0 * scale));
        } else {
            painter.add(Shape::line(curve.clone(), stroke));
        }
        let (Some(&tip), Some(&before)) = (curve.last(), curve.iter().rev().nth(2)) else {
            continue;
        };
        let along = (tip - before).normalized();
        let across = along.rot90();
        let size = (5.0 + 2.0 * strength) * scale.max(0.6);
        if bar {
            painter.line_segment([tip + across * size, tip - across * size], stroke);
        } else {
            painter.add(Shape::convex_polygon(
                vec![
                    tip + along * size,
                    tip + across * size * 0.7,
                    tip - across * size * 0.7,
                ],
                colour,
                Stroke::NONE,
            ));
        }
    }

    let font = |size: f32| egui::FontId::proportional(size * scale.max(0.7));
    let mut clicked = Clicked::Nothing;
    for (i, input) in layout.inputs.iter().enumerate() {
        let rect = boxed(i);
        let value = reads(input);
        painter.rect_filled(rect, rect.height() / 2.0, fill(value));
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            format!("{}\n{value:.2}", label(input)),
            font(12.5),
            ink(value),
        );
        ui.interact(rect, ui.id().with(("input", i)), Sense::hover())
            .on_hover_text(format!(
                "{}: read from outside the cell, {value:.2} where it sits now.",
                label(input)
            ));
    }
    for (n, name) in names.iter().enumerate() {
        let rect = boxed(inputs + n);
        let value = levels.get(n).copied().unwrap_or(0.0);
        painter.rect_filled(rect, 5.0 * scale, fill(value));
        if chosen == Some(n) {
            painter.rect_stroke(
                rect.expand(2.0),
                6.0 * scale,
                Stroke::new(2.0, Color32::WHITE),
                egui::StrokeKind::Outside,
            );
        }
        painter.text(
            rect.center(),
            egui::Align2::CENTER_CENTER,
            format!("{}\n{value:.2}", label(name)),
            font(13.0),
            ink(value),
        );
        if let Some(Some((badge, on))) = badges.get(n) {
            painter.text(
                rect.center_bottom() + Vec2::new(0.0, 2.0),
                egui::Align2::CENTER_TOP,
                badge,
                font(11.5),
                if *on {
                    Color32::from_rgb(255, 210, 120)
                } else {
                    Color32::from_gray(150)
                },
            );
        }
        let response = ui
            .interact(rect, ui.id().with(("node", n)), Sense::click())
            .on_hover_text(explain(n));
        if response.clicked() {
            clicked = Clicked::Node(n);
        }
    }
    clicked
}

/// The colour of a level, the ramp the tissue's protein view uses, so that
/// a box and the cells coloured by the same protein match.
fn fill(value: f64) -> Color32 {
    let [r, g, b, _] = paint::bytes(paint::ramp(value as f32));
    Color32::from_rgb(r, g, b)
}

/// Text that reads on that colour: dark on the bright end of the ramp.
fn ink(value: f64) -> Color32 {
    if value > 0.72 {
        Color32::from_gray(20)
    } else {
        Color32::from_gray(235)
    }
}

/// How a protein or an input is written for a reader rather than in a
/// file: the usual name where there is one, the file's name otherwise.
pub fn label(name: &str) -> String {
    match name {
        "oxygen" => "O₂",
        "egf" => "EGF",
        "crowding" => "crowding",
        "stress" => "stress",
        "waste" => "waste",
        "necrotic_factor" => "necrotic f.",
        "egfr" => "EGFR",
        "ras" => "RAS",
        "raf" => "RAF",
        "mek" => "MEK",
        "erk" => "ERK",
        "hif1" => "HIF-1α",
        "vegf" => "VEGF",
        "bcl2" => "BCL2",
        "bax" => "BAX",
        "mdm2" => "MDM2",
        "notch" => "NOTCH",
        "rb" => "RB",
        "e2f" => "E2F",
        "cyclin_d" => "cyclin D",
        "cyclin_e" => "cyclin E",
        "pdl1" => "PD-L1",
        other => other,
    }
    .to_owned()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn link(node: usize, source: &str) -> Link {
        Link {
            node,
            kind: "activated_by",
            source: Some(source.to_owned()),
            rate: 1.0,
            now: 0.0,
        }
    }

    #[test]
    fn a_relay_steps_one_column_further_and_a_loop_back_does_not_push_it() {
        // egf -> a -> b -> c, and c back onto b.
        let names: Vec<String> = ["a", "b", "c"].iter().map(|s| (*s).to_owned()).collect();
        let links = vec![link(0, "egf"), link(1, "a"), link(2, "b"), link(1, "c")];
        let layout = Layout::of(3, &names, &links);
        assert_eq!(layout.inputs, vec!["egf"]);
        assert_eq!(layout.columns, vec![0, 1, 2, 3]);
    }

    #[test]
    fn a_node_that_reads_nothing_sits_beside_the_inputs() {
        let names: Vec<String> = ["a", "b"].iter().map(|s| (*s).to_owned()).collect();
        let links = vec![link(1, "a")];
        let layout = Layout::of(2, &names, &links);
        assert!(layout.inputs.is_empty());
        assert_eq!(layout.columns, vec![1, 2]);
    }
}
