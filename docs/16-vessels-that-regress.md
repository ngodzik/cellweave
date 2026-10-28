# 16. Vessels that regress

In chapter 14 a sprout that never joined anything stayed for good, and a tumour whose sprouts had to join before they carried blood called in far more vessel than one fed at once. This chapter takes back the sprouts that carry no blood, as a real vessel network does, and asks whether that result was the model's or the tumour's.

## In the body

**A network is pruned.** A growing vessel network makes more connections than it keeps. **Pruning** removes segments, and **regression** is a segment's disappearance: its lumen, the channel blood flows through, collapses, and its endothelial cells (chapter 14) move into the vessels that stay. Korn and Augustin (Developmental Cell, 2015) review what is known of how it happens. Franco and colleagues (PLoS Biology, 2015) followed every endothelial cell of a remodelling network in the mouse retina, and found that "regression is mediated by dynamic and polarized migration of endothelial cells, representing anastomosis in reverse", and proposed that flow "leads to stabilization of high flow/shear segments and regression of adjacent low flow/shear segments". Lenard and colleagues (PLoS Biology, 2015) watched it live in the zebrafish embryo, where "blood vessel regression is a coordinated process of cell rearrangements involving lumen collapse and cell-cell contact resolution", and Chen and colleagues (PLoS Biology, 2012) found, in the zebrafish brain, that "changes to brain blood flow drive vessel pruning". A segment that carries little blood is the one that goes.

**In a tumour too.** Holash and colleagues (Science, 1999) found that some tumours begin by growing along the vessels already there, which is **co-option**, and that those vessels then regress, "leading to a secondarily avascular tumor and massive tumor cell loss", until "the remaining tumor is rescued by robust angiogenesis at the tumor margin". Vessels in a tumour are made and lost; what the tumour ends with is the balance of the two.

## In cellweave

**A sprout that carries nothing is taken back.** A configuration whose sprouts carry blood only once joined (chapter 14) may add:

```toml
[angiogenesis]
perfusion = "on_joining"
regress_after = 100   # steps a sprout may go without blood
```

A sprout that still carries no blood `regress_after` steps after it started regresses: every pixel of vessel it laid is given back to the tissue, which may move in. Every sprout grown from it goes with it, since blood reaching one of them would have reached the sprout it grew from, and without that sprout it is joined to nothing. A pixel it laid inside a disc a sprout still standing laid too stays, and belongs to that one. The vessel laid at the start never regresses, and neither does a sprout once it carries blood. The wall a regressed sprout left may sprout again, if the call is still there. The check comes after the tips have moved, so a sprout that joins another vessel on its last step stands.

Nothing is drawn, so a run with `regress_after` is the run without it until the first sprout regresses. Without it a blind end stays for good, as before. With `perfusion = "at_once"` every sprout carries blood from its first pixel and none could ever regress, so the configuration is refused.

**What a run says.** `measurements.csv` gains a `regressed` column, the sprouts taken back so far, which `Simulation::regressed` also says; a campaign writes it for each run, and a run's summary and the viewer's vessel line count them when sprouts can regress.

`examples/cord-angiogenic-regressing.toml` is `cord-angiogenic-joining.toml` with `regress_after = 100`.

```bash
cellweave --config examples/cord-angiogenic-joining.toml --versus examples/cord-angiogenic-regressing.toml --seeds 1-10 --mcs 300
```

### The predictions, and what happened

Written before anything was run with it, all at step 300 like chapter 14's comparison, seeds 1 to 10: (a) with `regress_after = 100`, fewer vessel pixels than without regression in at least 7 pairs of 10; (b) fewer living cells than without regression in at least 6, the blind ends kept without regression being joined later and then feeding; (c) fewer living cells than with blood at once in at least 6, taking away chapter 14's advantage of waiting for loops; (d) the median vessel at step 300 rising with the delay, from 50 steps to 100, to 200, to never; (e) at 100, at least 2 sprouts regressed by step 300, in the median.

| seed | living, at once | once joined | regressing after 100 | after 50 | vessel pixels, at once | once joined | after 100 | after 50 | regressed, after 100 | after 50 |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 47 | 88 | 88 | 63 | 185 | 1424 | 1424 | 466 | 0 | 6 |
| 2 | 96 | 105 | 78 | 49 | 508 | 1409 | 1547 | 171 | 5 | 2 |
| 3 | 59 | 87 | 65 | 36 | 240 | 1840 | 1413 | 222 | 7 | 7 |
| 4 | 64 | 49 | 49 | 57 | 260 | 216 | 495 | 502 | 1 | 4 |
| 5 | 45 | 54 | 49 | 49 | 144 | 201 | 481 | 493 | 1 | 3 |
| 6 | 65 | 58 | 57 | 58 | 360 | 343 | 655 | 621 | 2 | 7 |
| 7 | 125 | 68 | 65 | 37 | 874 | 1700 | 645 | 210 | 6 | 6 |
| 8 | 92 | 106 | 106 | 62 | 493 | 1363 | 1356 | 504 | 2 | 5 |
| 9 | 41 | 38 | 58 | 36 | 124 | 144 | 540 | 276 | 1 | 5 |
| 10 | 45 | 63 | 67 | 48 | 154 | 2106 | 573 | 226 | 6 | 5 |

(a) **Failed**: fewer vessel pixels in 4 pairs, more in 5, the same in 1. The median fell from 1386 pixels to 650 because the large networks shrank (seeds 3, 7 and 10), but the tumours that had started only one or two sprouts (seeds 4, 5, 6 and 9) ended with more vessel, not less. Without regression, a blind sprout whose tip had stopped stayed, and kept the wall around where it started from sprouting again, since a new sprout must start 16 pixels from where any other did (chapter 14). Taken back, it freed that wall, which sprouted again: seed 9 went from one sprout to four. So a blind end had been doing something nobody wrote: holding a patch of wall silent for good. (b) **Failed**: fewer living cells in 5 pairs, more in 2, the same in 3, 65 against 65.5 in the median, with a few more necrotic cells, 22.5 against 19. (c) **Failed**: more living cells than with blood at once in 6 pairs of 10 and fewer in 4, exactly as without regression, and fewer necrotic cells in 8. The tumour still called in more vessel than one fed at once, 2.6 times as much in the median where it had called 5.5 times. (d) **Failed as written**: 371 pixels at 50 steps and 650 at 100, but at 200 nothing had regressed by step 300, since no blind sprout was 200 steps old yet, and the ten runs were those without regression to the last bit. (e) **Held**: 2 regressed in the median, from none to 7.

At a delay of 50 steps the advantage is gone: fewer living cells than without regression in 8 pairs of 10, 49 in the median, and fewer than with blood at once in 7.

Then the same pairs to step 800, with two more predictions written first: (f) at 100 steps, more vessel than without regression in at least 6 pairs of 10, the walls regression frees going on sprouting where stopped blind ends had frozen the network; (g) at 200 steps, at least one sprout regressed in at least 8 runs of 10. A run stops when its tumour reaches the edge of the box, so the runs of a pair do not all end at the same step.

(f) **Failed, just**: more vessel in 5 pairs, less in 5. What happened is clearer counted by fate. Without regression, 4 tumours of 10 (seeds 4, 5, 6 and 9) never reached the edge: they had started one to three sprouts, all before step 300, and held 37 to 61 living cells at step 800. With it, all 10 reached the edge, those four between steps 506 and 633, after starting 33 to 37 sprouts, of which 7 to 13 regressed. The six that reached the edge either way did so a little later with regression in five, by 12 to 47 steps, and 4 steps sooner in one. Taking blind ends back slowed the tumours that were growing anyway, and released the ones a stopped blind end had frozen. (g) **Held**: at least one sprout regressed in 9 runs of 10. And the same four tumours were released: seeds 4, 5 and 9 took back a single sprout each and seed 6 two, and all four reached the edge between steps 547 and 619, where without regression they never did. One blind end taken back was enough.

So chapter 14's result was partly the model's. A tumour whose sprouts must join before they carry blood still calls in more vessel than one fed at once, and still has fewer dead, when a blind sprout survives 100 steps; when it survives 50, the tumour does worse than one fed at once. Whether waiting for loops helps a tumour or hurts it depends on how long a sprout that carries nothing is kept, a number this model does not know. And over a longer run regression did something no prediction expected: it kept tumours growing that would otherwise have stopped, because a blind end that stays also keeps the wall it left silent.

## Choices, and why

**No blood, rather than little.** Real pruning follows how much blood a segment carries. Here there is no flow (chapter 14), only whether a vessel carries blood or not, so the rule takes back the segments that carry none, the ones flow would remove first. A segment carrying blood is never pruned here, however redundant.

**A delay, not at once.** A sprout carries no blood while it grows, by definition, so taking back every dry sprout at once would let none of them ever join anything. The delay is the time a sprout has to find another vessel. 100 steps is about the cycle of a cord cell, and is a placeholder; the experiment above runs 50 and 200 too.

**The wall may sprout again.** The endothelium a sprout left still reads VEGF, vascular endothelial growth factor, the call of chapter 14. Forbidding the wall to sprout again would make every regression permanent, which nothing in the biology above says.

**What grew from it goes with it.** A sprout grown from a blind end is attached to nothing else, and carries no blood either. Leaving it standing would leave pieces of vessel floating in the tissue.

## Where in the code

`crates/core/src/types.rs`: the `regress_after` field of `Angiogenesis`. `crates/engine/src/vessels.rs`: `Sprout`, when it started and the discs it laid; `Tree::regress`, which takes a sprout back with what grew from it; `Vessels::regress_dry`, which finds the sprouts due. `crates/engine/src/cpm.rs`: `CpmLattice::clear_obstacle`. `crates/engine/src/simulation.rs`: where `regress_after` is checked, and `Simulation::regressed`. `crates/engine/src/measure.rs`: the `regressed` column. `src/campaign.rs`, `src/main.rs`, `viewer/src/panel.rs`: where the regressed sprouts are counted.

`examples/cord-angiogenic-regressing.toml`.

## Still a placeholder

The delay of 100 steps. Pruning by whether a vessel carries blood, not by how much: with no flow in the model, a loop that carries blood is never pruned, however long or redundant. A sprout regresses all at once, where a real one shrinks back as its cells move into the vessel they came from.

## References

- Chen Q, Jiang L, Li C, et al. [Haemodynamics-driven developmental pruning of brain vasculature in zebrafish](https://doi.org/10.1371/journal.pbio.1001374). *PLoS Biology* 10, e1001374 (2012). [Free full text](https://europepmc.org/article/PMC/PMC3419171).
- Franco CA, Jones ML, Bernabeu MO, et al. [Dynamic endothelial cell rearrangements drive developmental vessel regression](https://doi.org/10.1371/journal.pbio.1002125). *PLoS Biology* 13, e1002125 (2015). [Free full text](https://europepmc.org/article/PMC/PMC4401640).
- Holash J, Maisonpierre PC, Compton D, et al. [Vessel cooption, regression, and growth in tumors mediated by angiopoietins and VEGF](https://doi.org/10.1126/science.284.5422.1994). *Science* 284, 1994-1998 (1999).
- Korn C, Augustin HG. [Mechanisms of vessel pruning and regression](https://doi.org/10.1016/j.devcel.2015.06.004). *Developmental Cell* 34, 5-17 (2015).
- Lenard A, Daetwyler S, Betz C, et al. [Endothelial cell self-fusion during vascular pruning](https://doi.org/10.1371/journal.pbio.1002126). *PLoS Biology* 13, e1002126 (2015). [Free full text](https://europepmc.org/article/PMC/PMC4401649).
