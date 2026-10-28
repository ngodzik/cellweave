# 17. Cells pushed aside

A vessel growing through a tumour here takes the pixels of the cells in its way, and a cell left with none dies. This chapter counts how many die that way, lets the vessel push them aside instead, and asks whether the results of chapter 14 depended on it.

## In the body

**A sprout makes its way.** A sprout grows through the **extracellular matrix**, the mesh of protein fibres between cells, and its tip cell (chapter 14) digests the matrix in front of it with enzymes it secretes. Bauer and colleagues (PLoS Computational Biology, 2009), whose speeds chapter 14 uses, model "proteolytic matrix degradation by the tip cell", the matrix broken down by those enzymes, and how the density of its fibres sets how fast a sprout can go.

**The tissue makes room.** Models that build vessels out of cells answer what happens to the cells in a sprout's way by mechanics. In the three dimensional Cellular Potts model of Shirinifard and colleagues (PLoS One, 2009), vessels are made of cells of the lattice, which move only by displacing their neighbours as every cell does, and the pressure runs both ways: "the growing tumor pushes and stretches the vessels it contacts". In neither does a cell die for being in a sprout's way. In chapter 14 it did because a disc of new vessel was laid over whatever pixels were there, a shortcut of the model, not something seen in a tissue.

## In cellweave

**Counting the crushed.** `Simulation::crushed` counts the cells a growing vessel left with no pixel at all. A crushed cell is resorbed the step it is crushed, so it never shows among the necrotic cells at the end of a run: it is a living cell fewer, and one more gone. A run's summary says how many, a campaign writes a `crushed` column, and for a pair of arms with vessels it tallies which had more cells crushed, and which more vessel at the end.

**Pushing instead.** A configuration may add:

```toml
[angiogenesis]
in_the_way = "pushed"   # the default is "crushed"
```

Then, where a tip lays a disc of vessel, each cell pixel of the disc is shoved one pixel outward, the nearest to the centre first: it moves along the line from the disc's centre through it, and so does every cell pixel beyond it on that line, up to the first pixel of medium, which the last of them fills. The centre itself moves the way the tip is going. No cell changes size: a row of cells shifts by a pixel each, and the tissue grows into the medium around it, as it would grow into stroma. Only where the line meets another vessel or the edge of the box before any medium is the pixel taken, as before. `CpmLattice::push_obstacle_disc` does it.

Nothing is drawn, so pushing changes nothing until a disc first covers a pixel of a cell. Without `in_the_way`, cells are crushed, and every run of chapters 14 to 16 is the same as before.

`examples/cord-angiogenic-pushing.toml` is `cord-angiogenic.toml` with `in_the_way = "pushed"`.

```bash
cellweave --config examples/cord-angiogenic.toml --versus examples/cord-angiogenic-pushing.toml --seeds 1-10 --mcs 300
```

### The predictions, and what happened

Written before anything was run with it, at step 300, seeds 1 to 10: (a) at once, at most 3 cells crushed by step 300 in the median; (b) with blood only through loops, more cells crushed than at once in at least 7 pairs, there being five times the vessel; and with `in_the_way = "pushed"`: (c) no cell crushed in at least 8 runs of 10; (d) the living count barely changed, pushed over crushed between 0.9 and 1.1 in the median; (e) still more living cells than the plain cord of chapter 14 in at least 9 of 10; (f) with blood only through loops, still more vessel than at once in at least 7 pairs.

| seed | living, plain cord | crushed | pushed | vessel pixels, crushed | pushed | sprouts, crushed | pushed |
|---|---|---|---|---|---|---|---|
| 1 | 38 | 47 | 48 | 185 | 187 | 1 | 1 |
| 2 | 39 | 96 | 169 | 508 | 1391 | 5 | 12 |
| 3 | 39 | 59 | 161 | 240 | 1279 | 2 | 8 |
| 4 | 39 | 64 | 129 | 260 | 916 | 2 | 7 |
| 5 | 37 | 45 | 46 | 144 | 144 | 1 | 1 |
| 6 | 36 | 65 | 103 | 360 | 568 | 3 | 5 |
| 7 | 39 | 125 | 103 | 874 | 556 | 6 | 5 |
| 8 | 36 | 92 | 123 | 493 | 723 | 4 | 7 |
| 9 | 37 | 41 | 40 | 124 | 124 | 1 | 1 |
| 10 | 38 | 45 | 45 | 154 | 154 | 1 | 1 |

(a) **Held**, and more than held: not one cell was crushed, in any of the 40 runs, crushed or pushed, at once or through loops. A disc of vessel took pixels from many cells, and never a cell's last one. (b) **Failed**: none against none. (c) **Held**, for the same reason, and so says nothing.

(d) **Failed**, and by far: pushed, the tumour had more living cells in 7 pairs of 10, 103 against 61.5 in the median, and more than twice as many in seeds 2, 3 and 4. It started more sprouts, 5 against 2 in the median, and ended with more vessel in 6 pairs, as much in 3. Seed 3 shows how. The two runs are one run until step 128, when a disc first covers a pixel of a cell. Crushed, the only sprout's tip stops near step 175, and the vessel stays at 240 pixels to step 300; pushed, the tip goes on, and the vessel has reached 1066 pixels by step 275, with 128 cells born against 83. Why the tip stops when the cells in its way lose pixels was not tested. One reading is that a crushed cell gives up the pixels it secreted VEGF, vascular endothelial growth factor, from, right where the tip is, while a pushed one keeps them, a pixel further on; another, that a shrunken cell must grow back before it can divide.

(e) **Held**, in all 10, and by more: the tumour whose vessels pushed held 2.75 times the plain cord's living cells in the median, where the one whose vessels crushed held 1.58 times. (f) **Held**: pushed, blood only through loops still called in more vessel than blood at once in 8 pairs, 1325.5 pixels against 562 in the median, and still had fewer necrotic cells in 9, 20 against 41.5. But it had more living cells in only 3 pairs, 65.5 against 103 in the median, where crushed it had more in 6.

So the shortcut killed nobody, and still mattered. Chapter 14's first result, that vessels which grow feed a tumour, holds pushed, and is larger. That a tumour whose sprouts wait for loops calls in more vessel and has fewer dead holds too. That it also ends with more living cells does not: it was a result of the shortcut.

## Choices, and why

**Along a straight line.** Shifting the cells on one line by a pixel each is the smallest displacement that makes room without changing any cell's size. It leaves a cell dented where the line crossed it, and the lattice smooths the dent over the following steps, as it smooths any shape (chapter 1).

**The nearest pixels first.** A pixel near the centre pushes outward through the rest of the disc. Pushed later, it would meet pixels already turned to vessel, and be taken.

**Taken where nothing can move.** A vessel does not move, and nothing lies beyond the edge of the box. A cell between a new vessel and either has nowhere to go; it loses the pixel as before.

**Crushed by default.** Every run of the earlier chapters stays what it was, and a campaign can put the two side by side under the same seeds.

## Where in the code

`crates/core/src/types.rs`: `InTheWay`. `crates/engine/src/cpm.rs`: `CpmLattice::push_obstacle_disc`, and `shove`, which moves one line of pixels. `crates/engine/src/vessels.rs`: where a tip lays its disc, pushing or crushing. `crates/engine/src/simulation.rs`: `Simulation::crushed`. `src/campaign.rs`: the `crushed` column and the two tallies. `src/main.rs`: the summary line.

`examples/cord-angiogenic-pushing.toml`.

## Still a placeholder

Displacement along straight lines, one pixel per line, rather than the tissue relaxing as a whole. Vessels push cells but cells never push vessels, where in Shirinifard and colleagues' model the tumour deforms them. The matrix the sprout digests is not there at all. The patch of chapter 15 and the regression of chapter 16 were not run again with pushing.

## References

- Bauer AL, Jackson TL, Jiang Y. [Topography of extracellular matrix mediates vascular morphogenesis and migration speeds in angiogenesis](https://doi.org/10.1371/journal.pcbi.1000445). *PLoS Computational Biology* 5, e1000445 (2009). [Free full text](https://europepmc.org/article/PMC/PMC2709079).
- Shirinifard A, Gens JS, Zaitlen BL, et al. [3D multi-cell simulation of tumor growth and angiogenesis](https://doi.org/10.1371/journal.pone.0007190). *PLoS One* 4, e7190 (2009). [Free full text](https://europepmc.org/article/PMC/PMC2760204).
