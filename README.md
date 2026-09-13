# cellweave

**cellweave runs a molecular hypothesis inside every cell, and shows what it does to
the tissue.**

## What that means

A *molecular hypothesis* is a claim about biology: which proteins act on which, in
which direction, and how fast. For example, low oxygen raises HIF-1a, HIF-1a raises
BCL2, and BCL2 keeps the cell alive. Stated in words it sounds plausible. Written as
equations it becomes something you can run.

A *tissue* here is not a sheet of paper. It is a population of cells sitting next to
each other, pushing on each other, and sharing whatever diffuses between them.

So the program takes a hypothesis stated at the scale of molecules inside one cell,
gives a copy of it to every cell in the population, lets the cells consume and
release into the space around them, and reports what the whole population ends up
doing. The question it answers is not "is this protein important" but "if this
network is what is really inside each cell, what shape does the tissue take".

## Why

Going from data to a therapeutic target today is largely statistical association:
this gene correlates with that outcome. What sits between the two, and is hard to
get at, is the mechanism. A simulation is one way to state a mechanism precisely
enough that it can be wrong.

It also does something an experiment cannot. It runs the counterfactual: the same
population, the same random seed, one parameter changed, so the difference between
the two runs is attributable to that parameter and to nothing else. And it can be
searched, which means asking not only whether a hypothesis holds but which
interventions on it produce a tissue that behaves differently.

Nothing here is validated, and none of it substitutes for an experiment. The output
of a run is a hypothesis worth testing, never a result.

## How it works

Three scales, coupled, in one program:

- **The cells.** A Cellular Potts model on a lattice, where a cell is a set of
  pixels rather than a point. Cells have a shape, a surface, and neighbours they
  adhere to more or less strongly.
- **Inside each cell.** A small system of ordinary differential equations over
  protein concentrations, integrated per cell. Each cell carries its own state.
- **Around the cells.** Reaction diffusion fields for what moves through the space
  between cells, oxygen first.

Behaviour has to follow from that model rather than from hardcoded rules. A cell
grows, loosens its grip on its neighbours, or dies because its internal state says
so, not because a flag was set. A necrotic core, if one appears, has to appear
because the cells in the middle ran out of oxygen, not because the code drew a
circle.

## Status

Early-stage personal project. Not production software, and not scientifically
validated.

I am a software engineer, not a biologist. I started this to learn multiscale
modelling by implementing it. I use AI assistance to write code faster and to pick up
Rust idioms along the way, and I review and own everything here, mistakes included.

Mature, validated tools already exist for this kind of work. If you need results
today, use CompuCell3D, Morpheus, PhysiCell or PhysiBoSS instead. cellweave overlaps
with them by design, because reimplementing established methods is the point of the
exercise. What differs is the engineering: Rust, a single binary, reproducible TOML
configs.

Features will be added progressively, and this README will follow as they land.

## What works today

- Builds clean, `cargo clippy -D warnings` passes, 12 tests green
- 2D lattice with Monte Carlo spin flips and a volume constraint
- RK4 solver for a small per-cell protein network
- Explicit finite-difference diffusion solver with fixed boundaries
- TOML configuration, JSON snapshot output per saved step
- CLI entry point

## What does not work yet

Listed plainly, because a green test suite is not the same thing as correct physics.

- **The contact energy in the lattice is wrong**, so surface tension does not behave
  correctly. This is the first thing to fix.
- **The three scales are not coupled.** The binary runs the lattice and the protein
  network side by side with hardcoded inputs, and the diffusion solver is not wired in.
  Until that is done, the sentence at the top of this file describes the intent and
  not yet the program.
- **No parallelism.** `rayon` is declared as a dependency and unused.
- **No cell division, no cell death.** Cell count is static.
- **No stability guard on the diffusion solver.** A bad configuration diverges silently.
- **No calibration, no validation, no visualization.**

## Building and running

Requires a recent Rust toolchain (edition 2024, so 1.85 or newer).

```bash
cargo build --release
cargo test --workspace

# Run with built-in defaults: 200x200 lattice, one tumor cell
./target/release/cellweave --mcs 500

# Or from a config file
./target/release/cellweave --config sim.toml
```

Snapshots are written as `output/snapshot_NNNNNN.json`.

## Architecture

Separation between layers is enforced by the crate dependency graph, so the compiler
rejects inverted dependencies rather than relying on convention.

```
cellweave-core      domain types and traits, zero internal dependencies
      ^
cellweave-engine    lattice, ODE solver, diffusion
      ^
cellweave-io        TOML config, snapshot output
      ^
cellweave (bin)     CLI
```

There is no AI or MCP dependency inside cellweave. External tools connect through the
CLI. See [CLAUDE.md](CLAUDE.md) for coding rules.

## License

MIT. See [LICENSE](LICENSE).
