# 2. What happens inside a cell

How does a cell know there is food outside, and decide to grow?

## In the body

A cell does not decide to divide on its own. It waits for permission, and the permission is a molecule. A **growth factor** is a small protein released by other cells that means, roughly, "there is room and food here, grow". It binds a **receptor**, a protein sitting across the cell membrane with one end outside and one inside, and that binding starts a chain of relays inside the cell. Each relay is a protein that, once activated, activates the next. The chain is called a **signalling pathway**, and it ends in the nucleus, where the last relay switches genes on or off.

The pathway in cellweave is the one most often broken in cancer:

- **EGF**, epidermal growth factor, first found because it makes skin cells multiply. The input.
- **RAS**, a protein that works as an on/off switch, named after the rat sarcoma virus it was discovered in. When a receptor is bound, RAS flips on.
- **RAF** and **MEK**, two intermediate relays, kinases both. A **kinase** is an enzyme that activates another protein by attaching a phosphate group to it.
- **ERK**, extracellular signal-regulated kinase, the last relay. Active ERK enters the nucleus and switches on the genes for growth and division.

A mutation in RAS, found in roughly one cancer in three, leaves the switch stuck on. ERK stays active with no growth factor present, and the cell divides without being told. Several drugs target individual links of this chain.

Two more proteins matter here, both about surviving hardship rather than growing:

- **HIF-1α**, hypoxia-inducible factor 1 alpha. Under normal oxygen it is made and destroyed continuously, and there is little of it. As oxygen falls, its destruction slows and it accumulates. It is the cell's hypoxia response: it switches on genes for making new vessels, for surviving on less, and for stopping the cell cycle.
- **BCL2**, B-cell lymphoma 2, a survival protein. It holds off the cell's self-destruct machinery. When it collapses, the cell dies.
- **VEGF**, vascular endothelial growth factor, is what HIF-1α makes the cell secrete to call for new blood vessels. It is computed here and, for now, goes nowhere.

## In cellweave

Each cell carries five numbers between 0 and 1, the activated fraction of each protein: RAS, ERK, HIF-1α, VEGF, BCL2. And five equations say how each one changes over time given the others and what the cell sits in.

```mermaid
flowchart LR
    EGF(["EGF<br/>outside, an input"]) --> RAS["RAS"]
    RAS --> ERK["ERK"]
    O2(["oxygen<br/>outside, an input"]) -.->|"less oxygen,<br/>more HIF"| HIF["HIF-1α"]
    HIF --> VEGF["VEGF<br/>goes nowhere yet"]
    ERK --> BCL2["BCL2<br/>survival"]
    HIF --> BCL2
    O2 -.->|"below 15%:<br/>damage"| BCL2
    ERK ==>|"target volume"| G["the lattice<br/>chapter 4"]
    BCL2 ==>|"below 0.25: dead"| F["fate<br/>chapter 4"]
    HIF ==>|"above 0.45: stops cycling"| F
```

The equations, from `examples/networks/ras-erk.toml`, with each term named:

| protein | changes as | meaning |
|---|---|---|
| RAS | `0.5 egf (1 - ras) - 0.3 ras` | switched on by EGF, switches itself off at a steady rate |
| ERK | `0.8 ras (1 - erk) - 0.4 erk` | switched on by RAS, likewise |
| HIF-1α | `0.6 (1 - o2) (1 - hif1) - 0.5 hif1` | made faster the less oxygen there is, destroyed at a steady rate |
| VEGF | `0.7 hif1 - 0.3 vegf` | made in proportion to HIF-1α |
| BCL2 | `0.2 erk + 0.1 hif1 - 0.15 (bcl2 - 0.5) - 0.25 anoxia` | fed by ERK and HIF-1α, relaxes toward 0.5, and damaged below the anoxia threshold |

Every `(1 - x)` is saturation: a protein already fully active cannot be activated further. This is why ERK settles near 0.55 even when drowned in EGF.

The last term of BCL2 is the one that kills. `anoxia` is zero above 15% oxygen and rises to 1 at zero oxygen. Below the threshold, damage pulls survival down regardless of what ERK and HIF-1α feed in. That is deliberate: **necrosis is not a decision, it is a failure**. Below a critical oxygen level the cell cannot make enough energy, damage accumulates, and no survival programme rescues it. It is distinct from **apoptosis**, the orderly, programmed death a cell can choose. Only necrosis exists here so far.

### The network is a file, not code

Those rate constants used to be literals in the Rust source. They are now a
file, and that is the difference between a program that simulates one hypothesis
and a program that simulates whichever one you hand it. Stating a different
hypothesis is writing a different file.

Each line of the table above is a **term**, and terms come from a closed
catalogue, each named after what it means in biology rather than after its
algebra:

| term | what it computes | what it means |
|---|---|---|
| `activated_by` | `rate * source * (1 - self)` | one relay switching on the next, saturating because a protein cannot be more than fully active |
| `activated_by_absence_of` | `rate * (1 - source) * (1 - self)` | accumulating when something is scarce, which is HIF-1α as oxygen falls |
| `produced_by` | `rate * source` | made in proportion, without saturating, for something secreted |
| `decays` | `- rate * self` | first order turnover, which is what gives a node a steady state at all |
| `relaxes_to` | `- rate * (self - value)` | held near a set point rather than driven to zero |
| `damaged_below` | `- rate * max(0, (threshold - source) / threshold)` | energy failure, the one term nothing above cancels |

`activated_by` also takes a **Hill** response, `source^n / (k^n + source^n)` in
place of `source`, where `n` is how sharply it switches and `k` the level at
which it is half way. With a large `n` a graded signal becomes a decision.

A closed catalogue rather than free-form expressions is a deliberate limit. An
expression language would be more expressive and would let anyone write
something that parses and means nothing; a catalogue says exactly which
mechanisms the model knows about, and every one of them is documented above and
pinned by a test. A paper whose network does not fit adds an entry, which is a
change on the record rather than an open hole.

Two more things the file declares. **Roles**: nothing in a network says which of
its nodes answers the questions the simulation asks a cell, so the file says
which node is growth, which is survival, which is the hypoxia response. Naming
them rather than guessing from a node's name is what lets a network use its own
vocabulary. And **mechanics**: how the growth node becomes a target volume and a
stiffness for the lattice.

A file may also **extend** another and give only what differs:

```toml
extends = "ras-erk.toml"

[[nodes]]
name = "ras"
initial = 0.05
terms = [
  { kind = "relaxes_to", value = 1.0, rate = 0.5 },
  { kind = "decays", rate = 0.3 },
]
```

That whole file is the mutant network: RAS switched on by nothing instead of by
EGF, which is what an activating mutation in RAS does. Everything else comes
from the file it extends.

Which matters for a reason beyond brevity. A hypothesis should be the lines that
state it, and a copy of seventy identical lines hides the one that counts. It
also removes a way to be wrong: adjusting a rate in the healthy network without
adjusting the copy would leave two runs differing by two things, and the
comparison would mean nothing.

A node whose name appears in both **replaces** the inherited one whole, rather
than merging term by term. Replacing is the more predictable of the two: what a
node is, after the override, is written in one place. A new name is appended,
and roles and mechanics are inherited when left out. A file that extends itself
round a loop is refused, with the loop named.

Every name a term reads is resolved when the file is loaded, against the nodes
and against the two names the environment provides, `oxygen` and `egf`. A term
reading anything else is an error naming what it read and listing what is known,
rather than a silent zero in the middle of a run.

### How the equations are solved

Five equations that say "the rate of change of x is this function of x and the inputs" are an **ODE** system, ordinary differential equations: they involve time only, no space. They are advanced with **RK4**, the fourth order Runge-Kutta scheme, the standard way to integrate an ODE numerically with good accuracy. Each step of one unit of network time, the scheme evaluates the rates four times and blends them.

The network settles in about fifteen units of its time. The binary gives it eight units per lattice step, so it answers in about two. Why eight is chapter 3.

## Choices, and why

**Concentrations, not molecules.** Each protein is one number, a fraction between 0 and 1, as if the cell were a well-stirred bag. A real cell has thousands of copies of each protein in specific places. Averaging is the standard first model and it keeps five equations rather than five million.

**A reduced chain.** RAF and MEK are skipped: RAS activates ERK directly. They relay without deciding, and two more equations would add constants without adding behaviour. If a drug that targets MEK ever matters, they come back.

**ERK sets the target volume**, with the coefficient 1.5: `target = base (1 + 1.5 erk)`. It was 0.5 at first, and since ERK saturates near 0.55, no cell could ever reach the division threshold. The coefficient is chosen so a well stimulated cell aims at about twice its reference size. A model constant that silently makes a behaviour impossible is the kind of thing to state in a comment, not tune quietly.

**The damage term is separate from HIF-1α.** HIF-1α is pro-survival, and in this network hypoxia raises BCL2. That is real biology, and it means that without a damage term no cell could ever die of hypoxia. The damage term acts on BCL2 directly and is not offset by anything, because energy failure is not something a transcription programme can argue with.

## Where in the code

`crates/core/src/network.rs` is the hypothesis as written, and computes nothing:
`Term` is the catalogue above, `NodeSpec` a protein, `Roles` which node answers
which question, `MechanicsSpec` how growth reaches the lattice, `NetworkSpec` the
whole thing.

`crates/engine/src/network.rs` is it made runnable. `Network::from_spec` resolves
every name once, so the hot loop walks indices and never a string. `derivatives`
is the catalogue evaluated, `rk4` advances one step, `mechanics` turns the growth
node into a demand on the lattice, and `survival` and `reference_volume` are what
the binary reads to decide a fate.

`crates/io/src/config.rs` finds the file: `load_network` resolves it next to the
configuration, or falls back to `BUILT_IN_NETWORK`, which is this same RAS/ERK
file compiled in so that the built-in hypothesis and one you write go through
exactly the same path.

`crates/core/src/network.rs` also holds `merge_onto`, which is the rule above:
replace a node of the same name, append a new one, inherit what is left out.
`crates/io/src/config.rs` resolves the chain of files and refuses a loop.

`examples/networks/ras-erk.toml` is the network above, and
`examples/networks/ras-erk-mutant.toml` extends it with one node.

`crates/core/src/traits.rs` has `SignalingNetwork`, the contract any other network would have to meet. To simulate a different paper's network is to write another implementation, not to touch the engine.

## Still a placeholder

Every rate constant in the file is a plausible number, not a measured one. The literature has measured values for parts of this pathway, and calibrating against them is future work. What is fixed by the current numbers: ERK saturates near 0.55, the network settles in about fifteen units of its time, survival collapses between 5% and 6% of the bath's oxygen and is intact above 10%. Those are the facts the rest of the model was built on, and they are pinned by tests, so changing a constant that moves them will say so.
