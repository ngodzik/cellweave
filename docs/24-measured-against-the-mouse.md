# 24. Measured against the mouse

Chapter 12 set the sheet's clock by one measured rate and checked the shape of its clone sizes. This chapter sets the sheet's clones against the clones counted in the mouse oesophagus, one by one, at the times they were counted, and asks how far apart they are and what it takes to bring them together.

## In the body

**Clones counted, cell by cell.** Piedrafita and colleagues (Nature Communications, 2020) marked a few basal cells of the mouse oesophagus at random, each with a heritable colour, and counted their clones 10 days, 30 days, 84 days and 6 months later: for every clone, how many of its cells were still in the basal layer and how many had moved above it. They published every clone they counted, 1,222 of them from 16 mice, as a table under the article's licence, CC-BY 4.0, which lets anyone reuse it with credit.

| days after marking | clones counted | clones with a basal cell | mean basal cells in those |
|---|---|---|---|
| 10 | 300 | 250 | 2.07 |
| 30 | 315 | 281 | 4.37 |
| 84 | 302 | 295 | 6.83 |
| 182 | 305 | 298 | 11.18 |

**What the single progenitor model makes of them.** The model of chapter 12 fits these clones with a division rate of 2.9 per week, one division every 2.4 days, measured separately by the fading of a fluorescent protein; *r* = 0.1; and a share of basal cells that are progenitors, still dividing rather than committed to leave, of 0.65. In the authors' words: "The outcome of a given cell division is unpredictable but, on average, the likelihood of producing proliferating and differentiating cells is equal."

**Why *r* sets how fast clones grow.** A division whose daughters do opposite things, one dividing again and one leaving, adds a cell and then takes one away: the clone ends the size it was. Only a division whose daughters do the same thing changes it, by one either way. With *r* = 0.1 a fifth of divisions do, and a surviving clone grows by about *r* divided by the progenitor share, 0.15 cells, per division time. Over 84 days, about 35 divisions, that is five cells.

**Replacement next door.** Mesa and colleagues (Cell Stem Cell, 2018) followed every basal cell in large regions of the skin of living mice and found that "differentiation of a single stem cell is followed by division of a direct neighbor, but not vice versa". A cell leaving makes the room a neighbour divides into, which is the only coupling the sheet has (chapter 7).

## In cellweave

**Marking.** A configuration may add:

```toml
[tracing]
at = 1000    # the step every cell is marked at
every = 10   # every how many steps the marked clones are written down
```

At step `at`, before anything else happens, every cell in the layer is given a mark of its own, which its daughters inherit. A cell counts as in the layer as an experiment counts a basal cell: cycling, quiescent, or committed to differentiate and not yet gone. From then on, every `every` steps, `traced.csv` gains one row per mark some cell in the layer still carries: the step, the mark and how many cells carry it (`Simulation::traced_clones`, `TraceLog`). Marking draws nothing from the run's generator, so a traced run is the untraced run, step for step, which a test checks.

**The mouse's clones, fetched rather than kept.** The repository holds no one else's data (`CONTRIBUTING.md`). It holds a recipe, `datasets/piedrafita-2020.toml`: the paper, its licence, where the authors published the table (their Supplementary Data 10), and the fingerprint the file must have. `scripts/fetch-dataset.py` follows it: it downloads the file from the publisher, refuses it if the fingerprint differs, and writes one row per clone, the days since marking, the mouse, the clone, and its cells in the basal layer, above it, and terminally differentiated, to a table under `data/`, which git ignores. Whoever runs it fetches the clones for their own use, under the authors' licence.

**Comparing.** `scripts/compare-clones.py` takes runs and the mouse's table and, at each of the four times, sets the clones with a cell in the layer against the clones with a basal cell: their mean size, the shares larger than one, two and three times it, and the **Kolmogorov-Smirnov distance**, the largest gap between the two shares of clones at or below any size, 0 for samples from the same distribution and 1 for two that do not overlap. With `--division-days 2.414` each run's step stands for whatever makes its own division interval after marking 2.414 days, the mouse's, which is chapter 12's calibration done run by run: a setting that changes how fast the layer turns over is still compared after as many divisions as the mouse's. With `--progenitors 0.65` it also gives the share of the layer not committed to differentiate, averaged from marking on.

`examples/epithelium-traced.toml` is the neutral sheet of chapter 7, marked at step 1000. The tiled sheet starts with every cell the same size, and divides in waves for its first few hundred steps. The campaigns below ran it with NOTCH as chapter 7 set it; since this chapter, the sheet's default is the setting they found nearest the mouse (see the end of this section), so the same command now runs that one.

### The sheet as it was

Explored first, one run under the default seed gave surviving clones of 4.9, 12.0 and 36.6 cells at 10, 30 and 84 days: very nearly one plus the divisions elapsed, which is what clones do when every division and every departure is independent of every other, an effective *r* of one half. Predicted before the campaign, over seeds 1 to 10 pooled: (a) a larger mean than the mouse's at every time, by at least twice at 30, 84 and 182 days; (b) a larger ratio at 84 days than at 10; (c) a Kolmogorov-Smirnov distance of at least 0.3 at every time.

```bash
python3 scripts/fetch-dataset.py datasets/piedrafita-2020.toml
for s in $(seq 1 10); do
  ./target/release/cellweave --config examples/epithelium-traced.toml --seed $s --output runs/traced/$s
done
python3 scripts/compare-clones.py --data data/piedrafita-2020/lrig1-clones.csv \
  --at 1000 --division-days 2.414 --progenitors 0.65 runs/traced/*
```

| days | mean, mouse | mean, sheet | ratio | distance |
|---|---|---|---|---|
| 10 | 2.07 | 4.89 | 2.36 | 0.42 |
| 30 | 4.37 | 14.36 | 3.28 | 0.39 |
| 84 | 6.83 | 44.46 | 6.51 | 0.66 |
| 182 | 11.18 | 82.17 | 7.35 | 0.95 |

All three **held**. The share of the layer not committed was 0.88, against the mouse's 0.65.

The sheet's clones grow about six times as fast as the mouse's, per division. Chapter 12 found the sheet's sisters doing opposite things more often than chance, *r* = 0.22, and the clones do not show it: they grow as if *r* were a half. The sister count asks which of a daughter's next two events comes first; a daughter whose next event is leaving may leave a hundred steps after its birth, by which time its clone has had all that while to gain and lose cells elsewhere. A division leaves a clone's size alone only if the daughter that leaves does so soon after it is born, and nothing in the sheet makes it.

### Making replacement local

If the room a cell leaves were spread over the whole sheet by the volume constraint, whoever divided to fill it could be anywhere, and clones would mix as in a well stirred population. `scripts/voter-model.py` gives the other extreme, a **voter model**: every cell that leaves is replaced by a daughter of one of its six neighbours. Its clones hold 4.5, 10.2, 25.4 and 45.9 cells after 4, 12, 35 and 75 divisions per cell, the divisions of 10, 30, 84 and 182 days: replacement next door, alone, cannot bring the sheet near the mouse.

Predicted before trying it, over eight settings of ten seeds, temperature 10 or 3, `lambda_volume` 50 or 200, crowding as in chapter 7 or ten times stronger and half way at 75% of the reference size: (a) at 84 days a colder sheet gives smaller clones; (b) no setting gets them under twice the mouse's; (c) stronger crowding does not make them smaller.

| setting | divisions every | not committed | ratio at 10, 30, 84, 182 days |
|---|---|---|---|
| as in chapter 7 | 96 steps | 0.88 | 2.36, 3.28, 6.51, 7.35 |
| stiffer | 96 | 0.88 | 2.48, 2.91, 5.92, 7.03 |
| colder | 95 | 0.88 | 2.30, 2.70, 5.69, 6.70 |
| stronger crowding | 23 | 0.54 | 1.97, 2.29, 3.27, 3.26 |
| stronger crowding, stiffer | 22 | 0.54 | 2.01, 2.16, 2.95, 3.23 |
| stronger crowding, colder | 22 | 0.53 | 2.09, 2.08, 2.92, 3.04 |

(a) **Held**, barely. (b) **Held**. (c) **Failed**: crowding did more than anything else. A daughter is born at about 70% of its reference size, squeezed, and with crowding that strong it is often a newborn that leaves, often one of two sisters. At temperature 3 the stiffer and the softer sheet ran the same run, to the last bit: so cold, almost no flip that raises the energy is ever accepted, and which flips lower it depends on whether a cell is above or below its target, not on how stiff it is.

### Crowding

Predicted before a sweep of crowding's rate (0.3, 1, 3), its half way point (75% or 85% of the reference size, `k` 0.25 or 0.15) and NOTCH's baseline (as in chapter 7, or a fifth of it): (a) some setting under twice the mouse's clones at 84 days with between 0.45 and 0.8 of the layer not committed; (b) the clones smaller as crowding's rate rises; (c) the best setting within a distance of 0.2 at every time.

Ten settings of twelve emptied the layer: every cell committed and none was left dividing. Only a rate of 0.3, half way at 75%, held. With chapter 7's baseline it is the stronger crowding above. With a fifth of it, now the sheet's default, `examples/networks/epithelium.toml`, run 11,000 steps by `examples/epithelium-traced.toml`:

| days | mean, mouse | mean, sheet | ratio | distance |
|---|---|---|---|---|
| 10 | 2.07 | 2.77 | 1.34 | 0.14 |
| 30 | 4.37 | 6.34 | 1.45 | 0.15 |
| 84 | 6.83 | 18.30 | 2.68 | 0.34 |
| 182 | 11.18 | 55.11 | 4.93 | 0.65 |

with 0.87 of the layer not committed and a division every 96 to 108 steps across the seeds. All three **failed**.

For a month this sheet's clones are near the mouse's: a third larger, with distributions only 0.14 and 0.15 apart. Then they pull away. Making crowding, rather than the baseline, the reason a cell leaves made the newborns the ones that leave, and that slowed the clones; but a newborn squeezed enough to leave has a sister squeezed as much, and the two often leave together, which is a symmetric division as surely as two that both divide. The sisters' *r* in this sheet is 0.23, as in chapter 12's, while its clones grow half as fast: the sister count does not see what the clones see.

### What it says

The sheet couples a cell's departure to a neighbour's division through the room it leaves, as Mesa and colleagues saw, and its clones grow about six times as fast as the mouse's per division. The mouse's clones need about four divisions in five to end with one daughter staying and one gone soon after. Being squeezed at birth supplies some of that for a month, and not for six. What would supply the rest, a signal between sisters, an asymmetry at division, or a committed cell's room going to its own sister, the clones alone do not say; each makes a prediction about something else that could be measured, and the sheet can now be set against it.

### The default since

Of every setting tried, the last one's clones come nearest the mouse's at every time, and nothing it does is further from the mouse than the sheet as it was: the share of the layer not committed is the same, 0.87 against 0.88. It is now the sheet's NOTCH, in `examples/networks/epithelium.toml`, and the lesions of chapter 13 make NOTCH at 0.7 and 0.3 of its rate from both its sources. Its cells divide about every 104 steps rather than 97, so a step of the sheet is about 0.56 hours, which every sheet configuration now states. Chapters 7, 11, 12 and 13 ran with the setting before, and say so.

## Choices, and why

**Every cell marked.** The experiment marks a few cells so that two clones of one colour seldom touch. In the sheet every mark is distinct, so marking every cell gives about three hundred clones from one box, and no two can be mistaken for one.

**Marked once the sheet has settled.** At step 1000, past the waves of division a tiled sheet starts with, as an experiment marks an adult mouse.

**Counted as basal cells are.** A committed cell still in the layer counts, as a committed cell still in the basal layer does under a microscope; so does the share of the layer not committed, to set against the progenitor share.

**Calibrated run by run.** On the division rate after marking, so that settings which change how fast the layer turns over are compared after the same number of divisions.

**Surviving clones only, pooled.** The experiment's clones with a basal cell, pooled over its mice, against the runs' clones with a cell in the layer, pooled over seeds. Six months is read as 182 days.

## Where in the code

`crates/core/src/types.rs`: `Tracing`. `crates/engine/src/simulation.rs`: `Scenario::tracing`, where every cell is marked and a daughter inherits its mark, `Simulation::traced_clones`. `crates/io/src/config.rs`: the `[tracing]` table. `crates/io/src/output.rs`: `TraceLog`, `traced.csv`. `src/main.rs`: where the clones are written down. `crates/io/tests/epithelium.rs`: that marking changes nothing else.

`datasets/piedrafita-2020.toml`, `scripts/fetch-dataset.py`, `scripts/compare-clones.py`, `scripts/voter-model.py`. `examples/epithelium-traced.toml`, `examples/networks/epithelium.toml`, whose NOTCH is the setting this chapter found nearest the mouse.

## Still a placeholder

A box of about 300 cells, in which the six months' surviving clones hold a sixth of the layer or more. A committed cell leaves in a dozen to twenty steps, under a day, where the mouse's stays in the basal layer for about 1.3 days before moving up, on Piedrafita and colleagues' fitted rate of 5.4 a week. The settings tried are a grid, not a fit; nothing was done to find the best setting between them. And one experiment, in one tissue of one animal.

## References

- Piedrafita G, Kostiou V, Wabik A, et al. [A single-progenitor model as the unifying paradigm of epidermal and esophageal epithelial maintenance in mice](https://doi.org/10.1038/s41467-020-15258-0). *Nature Communications* 11, 1429 (2020). [Free full text](https://europepmc.org/article/PMC/PMC7080751).
- Mesa KR, Kawaguchi K, Cockburn K, et al. [Homeostatic epidermal stem cell self-renewal is driven by local differentiation](https://doi.org/10.1016/j.stem.2018.09.005). *Cell Stem Cell* 23, 677-686 (2018). [Free full text](https://europepmc.org/article/PMC/PMC6214709).
