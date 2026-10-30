# 26. A tissue with depth

Every tissue so far was one slice. A spheroid in a dish is a ball, and chapter 25 could only set the model's disc against measured balls through a conversion, a disc forming its core at $\sqrt{2/3}$ of the radius a ball does for the same oxygen. This chapter gives the lattice a depth: the same model, the same code, in a volume.

## In the body

**Tissue has three dimensions.** A tumour is a volume of cells pressed against each other on every side. In a packed tissue a cell touches a dozen or more neighbours, where in a section it seems to touch about six. Oxygen comes in from every direction a surface faces.

**A ball is fed better than a disc.** What starves the centre of a spheroid is the balance between what diffuses in through its surface and what the cells inside take up. A ball has more surface for its volume than a disc has edge for its area: $3/R$ against $2/R$. For the same uptake, the oxygen at the centre of a ball of radius $R$ falls short of the bath by $u R^2 / 6D$, of a disc by $u R^2 / 4D$, with $u$ the uptake and $D$ the diffusion coefficient. A core forms where that shortfall reaches what the cells can live on, so a ball forms its core at $\sqrt{3/2}$, about 1.22 times, the radius a disc does.

**Some tissues are flat.** The basal layer of the oesophagus of chapters 7, 11 and 24 is one cell thick, and cells grown on a dish lie in a layer. For them a slice is the geometry, not an approximation. A volume is for what has bulk: spheroids, tumours around a vessel.

## In cellweave

A configuration's `depth`, one when left out, is the number of slices. A plane is a volume one slice deep, and both are run by the same code.

```
    a plane: one slice             a volume: slices stacked

         . # .                       slice z-1      slice z      slice z+1
         # o #                         . . .         . # .         . . .
         . # .                         . # .         # o #         . # .
                                       . . .         . # .         . . .
    o has 4 neighbours #          o has 6: four in its slice, one before, one after
```

- **The lattice** stacks its slices one below the other, the voxel at `(x, y, z)` on row `z * height + y` of the grid, so that the first slice is exactly where a plane's pixels always were. A voxel has six neighbours: the four in its slice, up, down, left, right, then the one in the slice before and the one in the slice after.
- **A Monte Carlo step** makes one attempt per voxel, at voxels drawn at random. In a volume they are drawn only inside the box around everything that is not medium, one voxel wider on every side, as many as it has voxels: outside it a voxel and all its neighbours are medium and a flip changes nothing, so every voxel that can change is tried as often as before.
- **The first cell** is a ball of radius 3 voxels, 123 voxels, where in the plane it is a disc of radius 8, about 200 pixels. Why smaller is below.
- **A cell divides across a plane** through its centre of mass, whose normal is drawn uniformly over the sphere with two draws, where in the plane it divides along a line at an angle drawn with one.
- **A vessel seen end on** is a cylinder through the whole depth of the box, the disc of the plane in every slice. The tissue around it is a cord as long as the box, and the front and back faces cut through it, as the plane cuts through a cord.
- **The bath** reaches in from all six faces of the box, spreading through the faces of voxels, never their edges or corners, as in the plane it spreads through the sides of pixels and never their corners.
- **The field** is relaxed as in chapter 3, each voxel balanced against its six neighbours, in the red-black order, the colour of a voxel being the parity of $x + y + z$. The over-relaxation factor is the optimum for the box being swept, $2 / (1 + \sin(\pi / n))$ with $n$ its longest side, where the plane caps it at 1.95.
- **Radii** are those of a ball with the same volume, and around a vessel of a cylinder as long as the box, where in the plane they are those of a disc with the same area.
- **What runs in a volume**: a tumour in a bath or around a vessel, mutating or not, exposed to a lethal stress or not. A capillary lying in the plane, a sheet, vessels that grow, immune cells and lesions run only in the plane for now, and a configuration asking for one of them in a volume is refused with a message saying so.

The first test is the formula above. A ball of radius 10 voxels taking up 0.006 in a bath holding one settles with its centre at 0.47, against the formula's 0.50: the lattice's ball is a little larger than the ideal one, and starves a little more. A tumour grown from one cell in a box 120 voxels on a side, at an uptake of 0.0005 chosen to bring a core inside the box, had a radius of 28 voxels at step 200 and no dead cell; at step 250 a radius of 50, and a core. The formula puts the core at 48.

## Choices, and why

**One engine, not two.** A second, three-dimensional copy of the loop would be a second model, free to drift from the first without a test noticing, which is why the loop lives in the engine and not in the binary (chapter 1). The plane is the volume one slice deep.

**A plane is run exactly as it was.** Everything a volume adds comes after what the plane does, and only when there is a second slice: the third coordinate is drawn only in a volume, the two neighbours across slices come after the four in the slice, the third term of the side of a dividing plane is added only in a volume. So in the plane the same numbers are drawn, the same neighbours visited and the same sums taken in the same order. Every example in `examples/` was run to its end, the spheroid to step 300, before and after, and every file each wrote is the same, byte for byte.

**Coarser in a volume.** A cell sixteen pixels across in the plane would be some two thousand voxels across a ball, and a spheroid as large as the measured ones, some 350 µm in radius, tens of thousands of cells: a lattice of billions of voxels. At a first cell of radius 3, a cell is about seven voxels across, a voxel some 2.5 to 3 µm, and the measured spheroids fit in a box under three hundred on a side, about twenty million voxels. The cells are blockier; they are still cells, holding their size and shape against their neighbours.

**The time unit is the geometry's own.** In a volume the cells reach the size they divide at in fewer steps than in the plane. A cell of 123 voxels has most of them on its surface, where a flip can add one, and a cell of 200 pixels has a quarter; with room around it, a cell in a volume grows to its target within a few steps and divides at the next check for division, which comes every 20 steps (chapter 4). While the tissue is small the cells divide together, check after check: 2048 living cells at step 200, a power of two, in the 120 voxel box, whose volume doubled about every 20 steps where the plane's doubles every 50. The waves fade once room, not growth, decides who divides. A step is calibrated against a measured doubling time, as in chapter 25, separately for each geometry; `hours_per_step` belongs to a configuration for that reason (chapter 12).

**Only where the tissue is.** A spheroid starts as one cell in a box sized for its end, and for hundreds of steps nearly all of the box is empty medium. Three things now keep to the box around the tissue: the field's sweeps (chapter 3), the attempts of a Monte Carlo step in a volume, and the bath's flood, which holds everything outside the box around the cells without visiting it, since everything there is open and reaches a face of the box along a line that never enters it. In the plane the first and the last give the same result to the last bit, and the spheroid's first 120 steps on one core went from 65 seconds to 11.5. In a volume the attempts are drawn from a smaller box, which is a different sequence of draws from the same process.

**The factor that speeds the field up.** Successive over-relaxation settles fastest at a factor just under two, closer the larger the region (chapter 3). The plane caps it at 1.95, a little under the optimum for its boxes, and its runs were all made with that cap, so it stays. Capped there, a ball 250 voxels across settles about ten times slower than at its optimum; a volume takes the optimum for the box it sweeps, which still stays under two, where the floored iteration converges. On the 120 voxel box a step took 25 to 90 sweeps where it took 100 to 140.

**Dividing reads the lattice once per round.** Splitting a cell used to scan the whole lattice for its pixels, once per cell, which a volume of twenty million voxels and thousands of divisions cannot afford. The pixels of every cell ready to divide are now listed in one pass, which is the same lists, since splitting a cell relabels only its own pixels.

**A ceiling on living cells, per configuration.** The guard against runaway growth was 2000 living cells, enough for anything in the plane. A tumour in a volume needs tens of thousands, so the ceiling is now a configuration's `max_cells`, 2000 when left out.

**What is refused is refused loudly.** Vessels that grow, flow and immune cells are written against vessels laid in the plane; a sheet is a layer; a lesion is placed by its centre in the plane. Rather than run them wrongly in a volume, a configuration asking for one is refused, with a message.

## Where in the code

- `crates/engine/src/cpm.rs`: `CpmLattice::with_depth`, the `depth` field and the stacked `grid`; `von_neumann` returns `Neighbours`, six at most; `monte_carlo_step` and `occupied_span`, the box it draws from in a volume; `add_ball`; `pixels_of` and `divide_across`, which `divide` now calls; `voxels`, `occupant_at`, `position`, `dims3`.
- `crates/engine/src/pde.rs`: `ScalarField::with_depth`; `relax`, with the factor of a volume; `half_sweep` and `free_span` in three dimensions; `hold_bath_at`, from six faces, flooding only the box around what is closed; `hold_where`, `set_uptake_at`, `set_production_at`, `value_at`.
- `crates/engine/src/coupling.rs`: every function reads the lattice and writes the field by flat index, in a plane or a volume.
- `crates/engine/src/simulation.rs`: `Scenario::depth` and `Scenario::max_cells`, `DEFAULT_MAX_CELLS`; `FIRST_CELL_RADIUS` and `FIRST_CELL_RADIUS_3D`; `refuse_in_a_volume`; the plane a cell divides across in `divide_ready_cells`; `touches_the_edge` on six faces.
- `crates/engine/src/measure.rs`: `radius_of`, the radius of a disc, a ball or a cylinder; `depth_from_vessels`, slice by slice.
- `crates/io/src/config.rs`: `depth` and `max_cells`.

## Still a placeholder

- A first cell of radius 3 voxels, set by what a computer can hold rather than by a cell.
- The contact energies, λ and the temperature are the plane's. A voxel has six neighbours where a pixel has four, and nothing has been tuned for it.
- In a volume the cycle is the twenty steps between checks for division more than the time a cell takes to grow, and the first divisions come in waves.
- The voxel's size, about 2.8 µm, from the size of the melanoma spheroids' cells (chapter 27) and nothing else.

## References

- Browning AP, Sharp JA, Murphy RJ, Gunasingh G, Lawson B, Burrage K, Haass NK, Simpson MJ. [Quantitative analysis of tumour spheroid structure](https://doi.org/10.7554/eLife.73020). eLife 10, e73020 (2021). [Full text on Europe PMC](https://europepmc.org/article/PMC/PMC8741212).
- Graner F, Glazier JA. [Simulation of biological cell sorting using a two-dimensional extended Potts model](https://doi.org/10.1103/PhysRevLett.69.2013). Physical Review Letters 69, 2013 (1992).
