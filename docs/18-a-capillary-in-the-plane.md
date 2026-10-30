# 18. A capillary in the plane

The vessel of the cord is seen end on: a disc in the middle of the box, with blood running through it out of the plane. Nothing can flow through the vessels that grow from it, because they all start and end on the same disc. This chapter lays a capillary across the tissue instead, with an end where blood comes in and one where it leaves, and asks what the tumour does against it.

**Calibrated since.** Chapter 25 fixed the tumour's oxygen uptake on measured melanoma spheroids at 0.00002, where this chapter's runs used 0.0001. The mechanisms are the same; the shipped examples now run at the calibrated value, and rerunning them gives other numbers than those quoted here.

## In the body

**Vessels run through tissue.** A capillary is a tube, and tissue lies along it as well as around it. Some tumours begin by growing along the vessels already there, which is **co-option**, as Holash and colleagues (Science, 1999) found (chapter 16), before they call any new ones. Seen in a section cut across a vessel, the tissue is a sleeve around a disc, the cord of chapter 5; cut along it, a layer on each side of a line.

**Blood needs somewhere to go.** Blood flows from where its pressure is high to where it is low. A loop of new vessel that starts and ends on the same vessel, at the same pressure, carries nothing, however open it is. Owen and colleagues (Journal of Mathematical Biology, 2009) built a model of tissue whose vessels grow, carry blood and are pruned, and stress the role of "the elimination of vessels that do not sustain flow". For any of that to happen in cellweave, the vessels need an inlet and an outlet. That is what this chapter adds; the next one lets blood flow.

## In cellweave

**A capillary across the box.** A configuration may say:

```toml
[oxygen_source]
kind = "capillary"
radius = 3   # half its width: the capillary is 7 pixels wide
```

The capillary runs along a row a quarter of the way down the box, `capillary_row`, from its left side to its right, held at full oxygen like the cord's vessel. One cell is seeded against its lower wall, and the tissue grows into a layer along it, with three quarters of the box below to grow into. Cells cannot cross it, so the strip above stays stroma. `OxygenSource::is_vessel` says whether oxygen comes from a vessel of either kind, and everything that held a vessel's wall at full oxygen, the oxygen, the stress, the vessels that grow, holds the capillary's.

**A window onto a longer capillary.** The box is a window cut out of tissue that goes on beyond its left and right sides, the way the plane is a slice cut through a cord. So a run stops only when its tissue reaches the top or the bottom of the box, not its sides.

**How deep the layer is.** For a capillary, `viable_depth` is the living area divided by the number of columns living cells hold: the mean thickness of the living layer.

`examples/capillary.toml` is the layer; `examples/capillary-angiogenic.toml` lets vessels grow from it.

### Exploring first

The first runs were not predictions but a search for a layer comparable with the cord, and they found two things.

**A line feeds deeper than a point.** At the cord's uptake, 0.00003, the living layer was 97 pixels deep, nearly twice the cord's sleeve of 54. Around a disc the oxygen spreads over a ring that widens with distance; along a line it spreads over a strip that does not, so it reaches further. At 0.0001, three times the cord's uptake, the layer was 55 pixels deep, and the capillary's examples use that.

**At the cord's threshold, no wall sprouts.** At a VEGF threshold of 0.5 (VEGF, vascular endothelial growth factor, chapter 14), not one wall of the capillary sprouted in 800 steps, at any of three uptakes. Why the call arrives weaker was not tested. One difference is that the cord's starving cells come to surround its small disc on every side, where the layer's lie on one side of a long wall. At 0.3 the default seed started 60 sprouts and reached the bottom of the box by step 508; the angiogenic example uses 0.3. That the threshold is one absolute level whatever the geometry is the same placeholder chapter 14 points at.

### The predictions, and what happened

Written before the campaign, after those runs: (a) with vessels that grow, more living cells than the plain layer in at least 8 pairs of 10, at step 400; (b) with vessels that grow, the tissue reaching the bottom of the box by step 400 in at most 2; (c) the plain layer's living tissue lying within 48 to 59 pixels of the capillary in the median, within a tenth of the cord's 53.6 (chapter 14).

```bash
cellweave --config examples/capillary.toml --versus examples/capillary-angiogenic.toml --seeds 1-10 --mcs 400
```

| seed | living, plain | with vessels that grow | necrotic, plain | with vessels | sprouts | living tissue within, plain (px) |
|---|---|---|---|---|---|---|
| 1 | 66 | 200 | 46 | 0 | 33 | 54 |
| 2 | 63 | 167 | 37 | 0 | 33 | 55 |
| 3 | 61 | 152 | 47 | 0 | 28 | 54 |
| 4 | 71 | 163 | 39 | 0 | 33 | 54 |
| 5 | 80 | 171 | 47 | 0 | 33 | 55 |
| 6 | 73 | 152 | 35 | 0 | 30 | 54 |
| 7 | 72 | 141 | 44 | 0 | 28 | 52 |
| 8 | 61 | 134 | 45 | 0 | 26 | 53 |
| 9 | 69 | 140 | 45 | 0 | 34 | 54 |
| 10 | 55 | 177 | 41 | 0 | 31 | 53 |

(a) **Held**, in all 10: 157.5 living cells against 67.5 in the median, 2.3 times as many. (b) **Held**: no run reached the bottom by step 400. (c) **Held**: 54 pixels in the median, from 52 to 55.

What no prediction asked about stands out more. With vessels that grow and carry blood from their first pixel, not one cell was necrotic at step 400, in any run, where the plain layer had 35 to 47. About 30 sprouts had started in each, and the living tissue lay within 28 pixels of a vessel, half the plain layer's depth: vessels perfused at once fed every cell before any could starve. The next chapter makes blood wait for the pressure that drives it.

## Choices, and why

**A quarter of the way down.** The tissue cannot cross the capillary, so a capillary across the middle would leave half the box to stroma. A quarter of the way down leaves three quarters to the tumour.

**The sides do not stop a run.** Along a capillary that runs out of the box at both sides, tissue reaching a side has reached the edge of the window, not of the tissue. Stopping there would stop every run as soon as the layer had spread along the capillary, which it does first.

**Three times the uptake.** Comparing the layer with the cord asks first that both reach as deep. The uptake is the number that sets the depth, and the one already a placeholder.

**One capillary, straight.** A real capillary bed is a mesh. One straight capillary is the least that gives blood an inlet and an outlet.

## Where in the code

`crates/core/src/types.rs`: `OxygenSource::Capillary` and `OxygenSource::is_vessel`. `crates/engine/src/simulation.rs`: `capillary_row`, where the capillary is laid and the first cell seeded, and where touching the sides does not stop a run. `crates/engine/src/measure.rs`: the layer's `viable_depth`. `viewer/src/panel.rs`: the capillary in the run's description.

`examples/capillary.toml`, `examples/capillary-angiogenic.toml`.

## Still a placeholder

The capillary's width, row and straightness; the uptake, tripled to match the cord's depth; the threshold of 0.3. The box's sides are walls the tissue does not cross, where the tissue they stand for would go on.

## References

- Holash J, Maisonpierre PC, Compton D, et al. [Vessel cooption, regression, and growth in tumors mediated by angiopoietins and VEGF](https://doi.org/10.1126/science.284.5422.1994). *Science* 284, 1994-1998 (1999).
- Owen MR, Alarcón T, Maini PK, Byrne HM. [Angiogenesis and vascular remodelling in normal and cancerous tissues](https://doi.org/10.1007/s00285-008-0213-z). *Journal of Mathematical Biology* 58, 689-721 (2009).
