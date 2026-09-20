# cellweave coding rules

## Language

**All code, comments, doc-comments, variable names, error messages, and commit messages must be in English.**

## Prose and Markdown

**Never use em dashes or en dashes in Markdown files (README, CLAUDE.md, docs).** Use a comma, a colon, parentheses, or a new sentence instead. This applies to all prose in the repository, not just the README.

**Never hard-wrap prose at a fixed column.** One paragraph is one line, one bullet is one line. Soft wrapping is the job of the renderer and the editor, and a fixed column makes every later edit reflow a whole block for nothing. Fenced blocks, mermaid diagrams and table rows keep their own line structure.

## The course

`docs/` is a course that follows the software, chapter by chapter, from the biology up. A pull request that changes the model, a mechanism, a constant with a stated reason, or a scenario, is not finished until the chapter that describes it says the new thing. Every symbol the course names must exist in the code; check with grep before committing. Define every acronym and protein name at its first appearance in a chapter.

## Architecture

Separation of concerns is **enforced by the crate dependency graph**, not by convention. The compiler refuses inverted dependencies.

```
crates/core    ← no internal dependencies
    ↑
crates/engine  (depends on core)
    ↑
crates/io      (depends on core + engine)
    ↑
cellweave      (bin, depends on core + engine + io)
```

| Crate              | Single responsibility                           |
|--------------------|-------------------------------------------------|
| `cellweave-core`   | Domain types + traits, zero I/O, zero compute   |
| `cellweave-engine` | CPM, per-cell ODE, PDE diffusion                |
| `cellweave-io`     | TOML config, JSON snapshot output                |
| `cellweave` (bin)  | CLI entry point                                 |

No AI/MCP dependency inside cellweave. External tools connect through the CLI.

---

## Idiomatic Rust patterns

### 1. Traits for polymorphic behavior
Traits are the only abstraction mechanism. No inheritance, no classes.

```rust
// In core: what a signaling network MUST do
pub trait SignalingNetwork: Send + Sync {
    fn step(&mut self, dt: f64, inputs: &CellInputs) -> &ProteinState;
    fn protein_state(&self) -> &ProteinState;
}

// In engine: a concrete implementation
pub struct HallmarksNetwork { /* ... */ }
impl SignalingNetwork for HallmarksNetwork { /* ... */ }
```

Use `impl Trait` for static dispatch (hot paths). Use `dyn Trait` only for heterogeneous collections (mixed cell types).

### 2. Newtypes for domain primitives
Never bare `u32` or `f64` where a domain type makes sense. The compiler prevents confusing a `Concentration` with an `Energy`.

```rust
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd, serde::Serialize, serde::Deserialize)]
pub struct CellId(pub u32);

#[derive(Debug, Clone, Copy)]
pub struct Concentration(pub f64); // mol/L

#[derive(Debug, Clone, Copy)]
pub struct Energy(pub f64); // CPM energy units
```

### 3. Typestate for simulation lifecycle
The compiler makes it impossible to call `step()` before `initialize()`.

```rust
pub struct Simulation<S: SimState> {
    inner: SimInner,
    _state: PhantomData<S>,
}

pub struct Configured;
pub struct Running;
pub struct Done;

impl Simulation<Configured> {
    pub fn start(self) -> Simulation<Running> { /* ... */ }
}

impl Simulation<Running> {
    pub fn step(&mut self) -> SimSnapshot { /* ... */ }
    pub fn finish(self) -> Simulation<Done> { /* ... */ }
}
// step() does not exist on Simulation<Configured>: compile error
```

### 4. Builder for configuration
Complex config → Builder with validation at `build()`.

```rust
SimulationConfig::builder()
    .grid(200, 200, 1)
    .temperature(10.0)
    .add_cell_type(CellType::tumor().target_volume(50.0))
    .signaling(HallmarksNetwork::default())
    .diffusion_field("O2", DiffusionParams { d: 0.25, decay: 0.001 })
    .build()? // validation here, returns Result
```

### 5. Iterators over cell populations
No indexed loops over cells. Use iterators and `rayon`.

```rust
// Sequential
cells.iter().map(|c| c.protein_state())

// Parallel (rayon), drop-in replacement
cells.par_iter().map(|c| c.protein_state())
```

---

## Non-negotiable rules

### Error handling
- `thiserror` in library crates: explicit error types
- `anyhow` in binaries (`main`)
- Never `.unwrap()` or `.expect("...")` in library code
- Never `.unwrap()` on domain logic `None`: return `Option` or `Result`

### Unsafe
- `#![deny(unsafe_code)]` in every crate, the binary included
- There is no `unsafe` anywhere. Lifting the deny in one crate is a decision to argue for in the pull request that needs it, not a standing permission

### Parallelism
- `rayon` exclusively: no manual threads, no `std::sync::Mutex` in hot paths
- CPM uses the checkerboard scheme: even and odd pixels updated alternately, guaranteeing no data races without locks

### Visibility
- `pub` only for items that are part of the crate's public API
- `pub(crate)` to share between modules within a crate
- No `pub use *`: re-exports are explicit

### Testing
- Unit tests `#[cfg(test)]` in the relevant module
- Integration tests in `crate/tests/`
- `proptest` for numerical invariants (energy conservation, mass conservation)
- Each engine: regression test on a documented reference scenario

### Documentation
- Every `pub` item: `///` doc-comment required
- Math formulas: `/// Computes $\Delta H = \sum_\sigma J(\tau, \tau')$`
- `# Examples` for non-trivial functions
- Comments explain WHY, not what the name already says

---

## Biological model

### Three coupled scales (core of the project)

```
Extracellular (PDE)   O2, VEGF, EGF, MMP
      ↕ inputs / secretions
Subcellular (ODE)     RAS, ERK, p53, HIF-1α, BCL2, E-cad (per cell)
      ↕ behavioral parameters
Cellular (CPM)        mechanics, adhesion, volume, morphology
```

### Upward coupling (ODE → CPM)
Protein state drives mechanical parameters:
- High ERK → increased target volume (proliferation)
- Low E-cadherin → reduced cell-cell adhesion (invasion)
- BCL2/BAX ratio → apoptosis threshold

### Downward coupling (PDE → ODE)
Extracellular fields feed into per-cell ODEs:
- Local O2 → HIF-1α (hypoxia response)
- Local EGF → EGFR → RAS (proliferation signaling)

### Extensibility
The protein network is an `impl SignalingNetwork`. To simulate a different paper: new trait implementation, no engine modification.

---

## Commands

```bash
cargo build --workspace
cargo test --workspace
cargo clippy --workspace --all-targets -- -D warnings
cargo fmt --all -- --check
```
