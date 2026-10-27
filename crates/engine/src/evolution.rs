//! Heritable variation: a daughter that is not quite a copy of its parent.
//!
//! Without variation nothing can be selected, and a tumour is above all a
//! population under selection. Each division may, with a probability the
//! scenario sets, scale one rate of the daughter's network; the daughter then
//! founds a clone, and everything descended from it belongs to that clone
//! until the next mutation. Whether a clone spreads or vanishes is not decided
//! here: it follows from what the changed rate does to its cells' fates, in
//! the place where they live.

use cellweave_core::{CellInputs, CloneId, Concentration};

use crate::network::Mutation;
use crate::simulation::{Fate, SimState, Simulation};

/// A clone: where it came from, when, and what changed.
#[derive(Debug, Clone, PartialEq)]
pub struct CloneRecord {
    /// Its identifier, which is also its position in the list of clones.
    pub id: CloneId,
    /// The clone it arose in, or `None` for the founding clone.
    pub parent: Option<CloneId>,
    /// The step at which its first cell was born.
    pub founded_at: u64,
    /// The change that founded it, or `None` for the founding clone.
    pub mutation: Option<Mutation>,
    /// The mutation was drawn but, in this run, not made: the clone exists,
    /// with its number, its parent and its cells, and its network is its
    /// parent's. See [`crate::simulation::Scenario::suppress`].
    pub suppressed: bool,
}

impl CloneRecord {
    /// The clone of the first cell, which every other clone descends from.
    pub(crate) fn founder() -> Self {
        Self {
            id: CloneId(0),
            parent: None,
            founded_at: 0,
            mutation: None,
            suppressed: false,
        }
    }
}

/// A cell is resistant when, under a given level of the stress, its
/// network's chance of starting apoptosis is below this fraction of what the
/// unmutated network's would be.
pub const RESISTANT_BELOW: f64 = 0.5;

/// Network time a response is given to settle before it is read: well past
/// the fifteen units a network of this kind takes (chapter 2).
const RESPONSE_TIME: u32 = 60;

impl<S: SimState> Simulation<S> {
    /// The share of living cells that are resistant to the stress at `level`,
    /// `None` when the network has no apoptosis role or nothing lives.
    ///
    /// Each living cell's network is asked how it would respond with full
    /// oxygen, no crowding and the stress at `level`, and so is the unmutated
    /// network; a cell whose chance of starting apoptosis is below
    /// [`RESISTANT_BELOW`] of the unmutated one's counts as resistant. The
    /// question is put to the network on a copy, out of the tissue, so that
    /// a cell sheltered in a refuge is not counted as resistant for being out
    /// of reach.
    pub fn resistant_fraction(&self, level: f64) -> Option<f64> {
        let responses = self.relative_responses(level)?;
        let resistant = responses.iter().filter(|&&r| r < RESISTANT_BELOW).count();
        (!responses.is_empty()).then(|| resistant as f64 / responses.len() as f64)
    }

    /// Each living cell's chance of starting apoptosis under the stress at
    /// `level`, as a fraction of the unmutated network's: one for a cell that
    /// responds as the first cells did, below one for a cell that has become
    /// less sensitive, by however much. `None` when the network has no
    /// apoptosis role. Asked of each network on a copy, out of the tissue, as
    /// [`Simulation::resistant_fraction`] explains.
    pub fn relative_responses(&self, level: f64) -> Option<Vec<f64>> {
        let inputs = CellInputs {
            o2: Concentration(1.0),
            egf: Concentration(0.2),
            vegf: Concentration(0.0),
            crowding: Concentration(0.0),
            stress: Concentration(level),
        };
        let yardstick = self.pristine().apoptosis_under(&inputs, RESPONSE_TIME)?;
        let mut responses = Vec::new();
        for (id, record) in self.lattice().cell_records() {
            if record.volume.0 == 0
                || !matches!(self.fate(id), Some(Fate::Cycling | Fate::Quiescent))
            {
                continue;
            }
            let Some(chance) = self
                .network(id)
                .and_then(|n| n.apoptosis_under(&inputs, RESPONSE_TIME))
            else {
                continue;
            };
            responses.push(chance / yardstick);
        }
        Some(responses)
    }

    /// Living cells in each clone, indexed by clone identifier. A clone whose
    /// cells are all dead or gone counts zero; it is extinct.
    pub fn clone_sizes(&self) -> Vec<usize> {
        let mut sizes = vec![0; self.clones().len()];
        for (id, record) in self.lattice().cell_records() {
            let alive = matches!(self.fate(id), Some(Fate::Cycling | Fate::Quiescent));
            if alive
                && record.volume.0 > 0
                && let Some(size) = self.clone_of(id).and_then(|c| sizes.get_mut(c.0 as usize))
            {
                *size += 1;
            }
        }
        sizes
    }
}
