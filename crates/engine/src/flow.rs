//! Blood flow through a network of vessels, as current through resistors.
//!
//! A network is junctions joined by segments. Each segment lets through a flow
//! proportional to the difference of pressure between its ends, the factor
//! being its conductance: by Poiseuille's law, for a tube of radius $r$ and
//! length $L$, $Q = \frac{\pi r^4}{8 \mu L} \Delta p$. Some junctions are held at
//! a pressure, where blood comes in and where it leaves; at every other one,
//! what flows in flows out. That is a linear system in the free pressures,
//! solved here by elimination, which for the few hundred junctions a growing
//! network has is exact and fast.
//!
//! Only the ratios matter to who carries more blood than whom, so units are
//! dropped: a conductance is $r^4 / L$ with lengths in pixels.

/// One vessel segment between two junctions.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct Segment {
    /// The junction at one end.
    pub a: usize,
    /// The junction at the other.
    pub b: usize,
    /// How much it lets through per unit of pressure: $r^4 / L$.
    pub conductance: f64,
}

/// The conductance of a tube of `radius` and `length` pixels, up to the
/// constant every tube shares.
pub(crate) fn conductance(radius: f64, length: f64) -> f64 {
    radius.powi(4) / length.max(f64::MIN_POSITIVE)
}

/// The pressure at every junction, given the ones held fixed.
///
/// A junction joined to no held one, through any path, has no pressure to
/// take and reads NaN; every segment touching it carries nothing.
pub(crate) fn pressures(nodes: usize, segments: &[Segment], held: &[(usize, f64)]) -> Vec<f64> {
    let mut pressure = vec![f64::NAN; nodes];
    let mut fixed = vec![false; nodes];
    for &(node, p) in held {
        pressure[node] = p;
        fixed[node] = true;
    }
    let mut neighbours: Vec<Vec<(usize, f64)>> = vec![Vec::new(); nodes];
    for s in segments {
        if s.a != s.b && s.conductance > 0.0 {
            neighbours[s.a].push((s.b, s.conductance));
            neighbours[s.b].push((s.a, s.conductance));
        }
    }
    // Which free junctions a held one can be reached from: only those have
    // a pressure. Breadth first from every held junction, in order.
    let mut reached = fixed.clone();
    let mut queue: std::collections::VecDeque<usize> = held.iter().map(|&(node, _)| node).collect();
    while let Some(n) = queue.pop_front() {
        for &(m, _) in &neighbours[n] {
            if !reached[m] {
                reached[m] = true;
                queue.push_back(m);
            }
        }
    }
    let free: Vec<usize> = (0..nodes).filter(|&n| reached[n] && !fixed[n]).collect();
    let mut index = vec![usize::MAX; nodes];
    for (i, &n) in free.iter().enumerate() {
        index[n] = i;
    }
    // Kirchhoff at each free junction: the sum over its segments of
    // conductance times (its pressure less its neighbour's) is zero.
    let k = free.len();
    let mut matrix = vec![0.0; k * k];
    let mut rhs = vec![0.0; k];
    for (i, &n) in free.iter().enumerate() {
        for &(m, g) in &neighbours[n] {
            matrix[i * k + i] += g;
            if fixed[m] {
                rhs[i] += g * pressure[m];
            } else {
                matrix[i * k + index[m]] -= g;
            }
        }
    }
    for (i, p) in solve(k, &mut matrix, &mut rhs).into_iter().enumerate() {
        pressure[free[i]] = p;
    }
    pressure
}

/// The flow through each segment, from its end `a` to its end `b`: negative
/// when blood runs the other way, zero where a pressure is missing.
pub(crate) fn flows(segments: &[Segment], pressure: &[f64]) -> Vec<f64> {
    segments
        .iter()
        .map(|s| {
            let q = s.conductance * (pressure[s.a] - pressure[s.b]);
            if q.is_finite() { q } else { 0.0 }
        })
        .collect()
}

/// Solve `matrix · x = rhs`, `k` by `k`, row by row, by Gaussian elimination
/// with partial pivoting. The matrix here is a weighted graph Laplacian with
/// every component tied to a held pressure, so it is never singular.
fn solve(k: usize, matrix: &mut [f64], rhs: &mut [f64]) -> Vec<f64> {
    for col in 0..k {
        let pivot = (col..k)
            .max_by(|&i, &j| {
                matrix[i * k + col]
                    .abs()
                    .total_cmp(&matrix[j * k + col].abs())
            })
            .unwrap_or(col);
        if pivot != col {
            for c in 0..k {
                matrix.swap(col * k + c, pivot * k + c);
            }
            rhs.swap(col, pivot);
        }
        let diagonal = matrix[col * k + col];
        if diagonal == 0.0 {
            continue;
        }
        for row in col + 1..k {
            let factor = matrix[row * k + col] / diagonal;
            if factor == 0.0 {
                continue;
            }
            for c in col..k {
                matrix[row * k + c] -= factor * matrix[col * k + c];
            }
            rhs[row] -= factor * rhs[col];
        }
    }
    let mut x = vec![0.0; k];
    for row in (0..k).rev() {
        let mut sum = rhs[row];
        for c in row + 1..k {
            sum -= matrix[row * k + c] * x[c];
        }
        let diagonal = matrix[row * k + row];
        x[row] = if diagonal == 0.0 { 0.0 } else { sum / diagonal };
    }
    x
}

#[cfg(test)]
mod tests {
    use super::*;

    fn segment(a: usize, b: usize, conductance: f64) -> Segment {
        Segment { a, b, conductance }
    }

    #[test]
    fn two_tubes_in_series_share_one_flow_and_split_the_pressure() {
        // 0 -- 1 -- 2, the second tube half as conductive as the first.
        let segments = [segment(0, 1, 2.0), segment(1, 2, 1.0)];
        let p = pressures(3, &segments, &[(0, 1.0), (2, 0.0)]);
        assert!((p[1] - 2.0 / 3.0).abs() < 1e-12, "{p:?}");
        let q = flows(&segments, &p);
        assert!((q[0] - q[1]).abs() < 1e-12, "what goes in comes out: {q:?}");
        assert!((q[0] - 2.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn a_bypass_takes_its_share_and_a_blind_end_takes_none() {
        // A capillary 0 -- 1 -- 2 -- 3 held at 1 and 0, a bypass from 1 to 2,
        // and a blind branch 1 -- 4 -- 5 that goes nowhere.
        let segments = [
            segment(0, 1, 1.0),
            segment(1, 2, 1.0),
            segment(2, 3, 1.0),
            segment(1, 2, 1.0),
            segment(1, 4, 1.0),
            segment(4, 5, 1.0),
        ];
        let p = pressures(6, &segments, &[(0, 1.0), (3, 0.0)]);
        let q = flows(&segments, &p);
        assert!((q[1] - q[3]).abs() < 1e-12, "two equal paths, equal flows");
        assert!((q[0] - (q[1] + q[3])).abs() < 1e-12, "and they add up");
        assert!(
            q[4].abs() < 1e-12 && q[5].abs() < 1e-12,
            "a blind end: {q:?}"
        );
        assert!((p[5] - p[1]).abs() < 1e-12, "at the pressure it hangs from");
    }

    #[test]
    fn a_piece_joined_to_nothing_held_has_no_pressure_and_no_flow() {
        let segments = [segment(0, 1, 1.0), segment(2, 3, 1.0)];
        let p = pressures(4, &segments, &[(0, 1.0), (1, 0.0)]);
        assert!(p[2].is_nan() && p[3].is_nan());
        assert_eq!(flows(&segments, &p)[1], 0.0);
    }

    #[test]
    fn conductance_follows_the_fourth_power_of_the_radius() {
        assert_eq!(conductance(2.0, 1.0), 16.0 * conductance(1.0, 1.0));
        assert_eq!(conductance(1.0, 4.0), conductance(1.0, 1.0) / 4.0);
    }
}
