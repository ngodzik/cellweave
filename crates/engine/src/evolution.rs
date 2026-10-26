//! Heritable variation: a daughter that is not quite a copy of its parent.
//!
//! Without variation nothing can be selected, and a tumour is above all a
//! population under selection. Each division may, with a probability the
//! scenario sets, scale one rate of the daughter's network; the daughter then
//! founds a clone, and everything descended from it belongs to that clone
//! until the next mutation. Whether a clone spreads or vanishes is not decided
//! here: it follows from what the changed rate does to its cells' fates, in
//! the place where they live.

use cellweave_core::CloneId;

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
}

impl CloneRecord {
    /// The clone of the first cell, which every other clone descends from.
    pub(crate) fn founder() -> Self {
        Self {
            id: CloneId(0),
            parent: None,
            founded_at: 0,
            mutation: None,
        }
    }
}

impl<S: SimState> Simulation<S> {
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
