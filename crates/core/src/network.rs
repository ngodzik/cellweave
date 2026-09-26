//! How a protein network is declared, before anything runs it.
//!
//! A network is a list of nodes and, for each, the terms that make its value
//! change. Nothing here computes: this is the shape a hypothesis takes on disk,
//! so that stating a different hypothesis is writing a different file rather
//! than editing and rebuilding the program.
//!
//! The terms are a closed catalogue, each named after what it means in biology
//! rather than after its algebra. That is deliberate. A free-form expression
//! would be more expressive and would let anyone write something that parses
//! and means nothing; a catalogue says exactly which mechanisms the model knows
//! about, and every one of them is documented and pinned by a test. A paper
//! whose network does not fit adds an entry here, which is a change on the
//! record rather than an open hole.

use serde::{Deserialize, Serialize};

/// What the environment offers a term as a source, besides the network's own
/// nodes.
///
/// These names are the engine's, not the configuration's: they exist because
/// something outside the cell computes them. More appear as the coupling grows.
pub const INPUT_NAMES: [&str; 2] = ["oxygen", "egf"];

/// Cooperative binding: `source^n / (k^n + source^n)` stands in for `source`.
///
/// `n` is how sharply the response switches, `k` the level at which it is half
/// way. With `n = 1` this is ordinary Michaelis-Menten saturation; large `n`
/// approaches a step at `k`, which is how a graded signal becomes a decision.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Hill {
    /// Cooperativity exponent.
    pub n: f64,
    /// Half-maximal level of the source.
    pub k: f64,
}

/// One term in the rate of change of one node.
///
/// Every term reads at most one source and writes one node. `self` below means
/// the node the term belongs to.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Term {
    /// `rate * source * (1 - self)`. The source switches this node on, and the
    /// node saturates because it cannot be more than fully active. The usual
    /// shape of one relay activating the next.
    ActivatedBy {
        /// Node or input that drives it.
        source: String,
        /// How fast, per unit time.
        rate: f64,
        /// Cooperative response instead of a proportional one.
        #[serde(default)]
        hill: Option<Hill>,
    },
    /// `rate * (1 - source) * (1 - self)`. The node accumulates when the source
    /// is scarce: HIF-1α as oxygen falls, whose destruction slows rather than
    /// whose production rises.
    ActivatedByAbsenceOf {
        /// Node or input whose scarcity drives it.
        source: String,
        /// How fast, per unit time.
        rate: f64,
    },
    /// `rate * source`. Made in proportion to the source and not saturating,
    /// for something secreted rather than switched on.
    ProducedBy {
        /// Node or input it is made from.
        source: String,
        /// How fast, per unit time.
        rate: f64,
    },
    /// `- rate * self`. First order turnover: a fixed fraction lost per unit
    /// time, which is what gives a node a steady state at all.
    Decays {
        /// Fraction lost per unit time.
        rate: f64,
    },
    /// `- rate * (self - value)`. Held near a set point rather than driven to
    /// zero, for something the cell maintains.
    RelaxesTo {
        /// The level it is held near.
        value: f64,
        /// How strongly.
        rate: f64,
    },
    /// `- rate * max(0, (threshold - source) / threshold)`. Damage below a
    /// critical level of the source, offset by nothing.
    ///
    /// This is how necrosis differs from every other term: energy failure is
    /// not a programme a transcription factor can argue with, so no amount of
    /// pro-survival signalling cancels it.
    DamagedBelow {
        /// Node or input whose shortage does the damage.
        source: String,
        /// Level below which damage begins.
        threshold: f64,
        /// How fast, at zero source.
        rate: f64,
    },
}

impl Term {
    /// The source this term reads, if it reads one.
    pub fn source(&self) -> Option<&str> {
        match self {
            Self::ActivatedBy { source, .. }
            | Self::ActivatedByAbsenceOf { source, .. }
            | Self::ProducedBy { source, .. }
            | Self::DamagedBelow { source, .. } => Some(source),
            Self::Decays { .. } | Self::RelaxesTo { .. } => None,
        }
    }
}

/// One protein, with where it starts and what moves it.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NodeSpec {
    /// How terms and roles refer to it.
    pub name: String,
    /// Its value before anything runs, between 0 and 1.
    #[serde(default)]
    pub initial: f64,
    /// Everything that makes it change.
    #[serde(default)]
    pub terms: Vec<Term>,
}

/// Which node plays which part, for the code outside the network.
///
/// The simulation asks a cell three questions: how large it is trying to be,
/// whether it is still alive, and whether it is short enough of oxygen to have
/// stopped cycling. Nothing in a network says which of its nodes answers those,
/// so the configuration does. Naming them rather than guessing from a node's
/// name is what lets a network use its own vocabulary.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Roles {
    /// Node whose level sets the target volume.
    pub growth: String,
    /// Node whose collapse means the cell is dead.
    pub survival: String,
    /// Node whose height means the cell has stopped cycling.
    pub hypoxia_response: String,
}

/// How the growth node becomes a mechanical demand on the lattice.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct MechanicsSpec {
    /// Target volume is `reference * (1 + coefficient * growth)`.
    ///
    /// Worth a moment: if the growth node saturates below `1 / coefficient`
    /// times the division threshold, no cell can ever divide. The first value
    /// tried in this project did exactly that.
    pub growth_coefficient: f64,
    /// Stiffness of the volume constraint handed to the lattice.
    pub lambda_volume: f64,
    /// Contact energy against medium handed to the lattice.
    pub j_medium: f64,
    /// Contact energy against a cell of the same kind, as
    /// `j_self_base + j_self_from_growth * (1 - growth)`.
    pub j_self_base: f64,
    /// How much a low growth node stiffens adhesion.
    pub j_self_from_growth: f64,
}

/// A whole hypothesis: the nodes, what each node is for, and how the network
/// reaches the lattice.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct NetworkSpec {
    /// The proteins, in the order they appear in a snapshot.
    pub nodes: Vec<NodeSpec>,
    /// Which node answers which question.
    pub roles: Roles,
    /// How the growth node becomes mechanics.
    pub mechanics: MechanicsSpec,
}
