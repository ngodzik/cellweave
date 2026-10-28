//! What two sisters did next.
//!
//! A division makes two cells, and each of them goes on either to divide
//! again or to leave. In a renewing epithelium the two are not independent:
//! lineage tracing in the mouse oesophagus and skin finds sisters more often
//! doing opposite things than chance would have them (Piedrafita and
//! colleagues, 2020). Nothing in the model says what a sister does, so whether
//! its sisters coordinate is something to measure, not something to set.

/// How one life, from the division that began it to the cell's next event,
/// ended.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Next {
    /// The cell divided again.
    Divided,
    /// The cell left: it differentiated, or died.
    Left,
}

/// Pairs of sisters, counted by what both did next. Only pairs whose two
/// next events have both happened are counted.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SisterFates {
    /// Both divided again.
    pub both_divided: u64,
    /// One divided and the other left.
    pub one_each: u64,
    /// Both left.
    pub both_left: u64,
}

impl SisterFates {
    /// Pairs counted.
    pub fn pairs(&self) -> u64 {
        self.both_divided + self.one_each + self.both_left
    }

    /// The share of divisions with a symmetric outcome on each side, the `r`
    /// of the single progenitor model: $r = (n_{PP} + n_{DD}) / 2n$.
    ///
    /// A quarter when sisters are independent and each is as likely to
    /// divide as to leave; below a quarter when sisters tend to do opposite
    /// things, above when they tend to do the same. `None` before any pair
    /// is complete.
    pub fn r(&self) -> Option<f64> {
        let pairs = self.pairs();
        (pairs > 0).then(|| (self.both_divided + self.both_left) as f64 / (2 * pairs) as f64)
    }
}

/// Every life a division began, its sister's, and how it ended.
#[derive(Debug, Clone, Default)]
pub(crate) struct Lives {
    /// For each life, the index of its sister's life, and how it ended once
    /// it has.
    lives: Vec<(usize, Option<Next>)>,
    /// The life the cell in each slot is living, `None` for a cell no
    /// division made, which has no sister.
    current: Vec<Option<usize>>,
}

impl Lives {
    /// No division yet, and `cells` cells that no division made.
    pub(crate) fn new(cells: usize) -> Self {
        Self {
            lives: Vec::new(),
            current: vec![None; cells],
        }
    }

    /// The cell in slot `parent` divided, and its daughter took slot
    /// `daughter`: the parent's life ends, and two sister lives begin.
    pub(crate) fn divided(&mut self, parent: usize, daughter: usize) {
        self.end(parent, Next::Divided);
        let first = self.lives.len();
        self.lives.push((first + 1, None));
        self.lives.push((first, None));
        if self.current.len() <= daughter {
            self.current.resize(daughter + 1, None);
        }
        self.current[parent] = Some(first);
        self.current[daughter] = Some(first + 1);
    }

    /// The cell in slot `slot` left the layer or died.
    pub(crate) fn left(&mut self, slot: usize) {
        self.end(slot, Next::Left);
    }

    fn end(&mut self, slot: usize, next: Next) {
        if let Some(life) = self.current.get_mut(slot).and_then(Option::take)
            && let Some(entry) = self.lives.get_mut(life)
        {
            entry.1 = Some(next);
        }
    }

    /// Every complete pair, counted once.
    pub(crate) fn fates(&self) -> SisterFates {
        let mut fates = SisterFates::default();
        for (i, &(sister, next)) in self.lives.iter().enumerate() {
            if sister < i {
                continue;
            }
            match (next, self.lives.get(sister).and_then(|s| s.1)) {
                (Some(Next::Divided), Some(Next::Divided)) => fates.both_divided += 1,
                (Some(Next::Left), Some(Next::Left)) => fates.both_left += 1,
                (Some(_), Some(_)) => fates.one_each += 1,
                _ => {}
            }
        }
        fates
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_pair_is_counted_once_both_its_lives_have_ended() {
        let mut lives = Lives::new(1);
        lives.divided(0, 1);
        assert_eq!(lives.fates().pairs(), 0, "neither has done anything yet");
        lives.left(1);
        assert_eq!(lives.fates().pairs(), 0, "one of two");
        lives.divided(0, 2);
        assert_eq!(
            lives.fates(),
            SisterFates {
                one_each: 1,
                ..SisterFates::default()
            }
        );
        // Slot 0's new sister is slot 2; both leave.
        lives.left(0);
        lives.left(2);
        assert_eq!(lives.fates().both_left, 1);
        assert_eq!(lives.fates().r(), Some(0.25));
    }

    #[test]
    fn a_cell_no_division_made_has_no_pair() {
        let mut lives = Lives::new(3);
        lives.left(2);
        lives.left(1);
        assert_eq!(lives.fates().pairs(), 0);
        assert_eq!(lives.fates().r(), None);
    }
}
