# cellweave

A multiscale cancer cell simulator written in Rust.

## Goal

Simulate how tumor cells behave by coupling three scales in a single program: the
shape and movement of the cells, the protein signaling inside each of them, and the
molecules diffusing around them.

The point is that behavior should follow from the model rather than from hardcoded
rules. A cell should grow, loosen its grip on its neighbors, or die because its
internal state says so, not because a flag was set.

## Status

Early-stage personal project. Not production software, and not scientifically
validated.

I am a software engineer, not a biologist. I started this to learn multiscale
modeling by implementing it. I use AI assistance to write code faster and to pick up
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
