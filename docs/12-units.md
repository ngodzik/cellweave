# 12. Units

How long a step is, measured against a mouse; two laws every renewing epithelium obeys, checked on the sheet; and one coordination between sister cells that nobody wrote, found in it.

## In the body

How a basal cell of the mouse oesophagus behaves has been worked out by two kinds of experiment, brought together by Piedrafita and colleagues (Nature Communications, 2020).

The first measures how fast cells divide. A fluorescent protein that sticks to DNA (chromatin) is made in every cell for a while, then switched off; each division halves what a cell carries, so how fast the glow fades says how often cells divide. In the mouse oesophagus, every basal cell divides at about the same rate: on average every 2.4 days. There is no slow, separate population of stem cells.

The second is lineage tracing (chapter 7), at low density: a few scattered basal cells marked, their clones counted months later. It shows the **single progenitor model** at work. Every dividing basal cell is alike, and a division makes two cells that will divide (with probability *r*), two that will differentiate (also *r*), or one of each (the rest). On average one of each, so the layer keeps its size, and which clone wins is chance. Two laws follow, and the mouse data obey both. Clones are lost steadily, so the survivors' mean size grows in proportion to time. And once enough time has passed, the sizes of the surviving clones, measured in multiples of their mean, settle into one fixed shape: the share of clones larger than *x* times the mean is $e^{-x}$, about 37% beyond the mean, 14% beyond twice it, 5% beyond three times.

The third finding is about sisters. If the two cells a division makes chose their fates independently, both would divide or both would leave at least half the time, and *r* would be a quarter or more. In the mouse oesophagus and skin, *r* is below a quarter: sisters do opposite things more often than chance. The authors put it down to local coordination between a cell leaving the layer and a neighbour dividing.

## In cellweave

Three things the run now reports, and one it can be told.

**How often a living cell divides.** The summary gives, over the second half of the run, the living averaged over that window times its length, over the births in it. In the neutral sheet of `examples/epithelium.toml` it is once every 94 to 102 steps across seeds 1 to 10, about 97.

**What sisters do next.** `Simulation::sister_fates` follows every division: the two cells it makes, and each one's next event, dividing again or leaving. A pair counts once both have happened. The summary gives the shares, and *r* (`SisterFates::r`), half the share of pairs that did the same thing.

**How the surviving clones spread.** For a sheet, the summary gives the share of surviving clones larger than one, two and three times their mean size, beside the exponential.

**Hours per step.** A configuration may state `hours_per_step`. Nothing in the simulation reads it; it converts what a run reports, and the summary says how long the run stood for. The sheet's configurations state 0.6: about 97 steps for a division here, 2.4 days for one in the mouse, so a step is about 36 minutes, the 3000 steps of `examples/epithelium.toml` about 75 days, and the 6000 of chapters 7 and 11 about 150 days of a mouse's life.

Space needs no unit beyond the cell. Clone sizes are counted in basal cells, here as in the papers, and the box holds 20 cells by 20.

### The predictions, and what happened

Written before any code, on the neutral sheet:

| prediction | result |
|---|---|
| from step 300 to 3000, clones still alive times the step stays within 20% of its value at step 300 (the default seed) | **held**: 29,400 at step 300, then between 26,000 and 34,200 |
| pooled over seeds 1 to 10 at step 1500, the shares larger than 1, 2 and 3 times the mean are within 0.07 of 0.37, 0.14 and 0.05 | **held**: 205 surviving clones; 0.34, 0.18 and 0.07 |
| sisters' fates anti-correlated: *r* below a quarter | **held**: 0.218 to 0.234 across the ten seeds, from about 4,200 pairs each by step 1500 |

The clone counts behind the first line, under the default seed:

| step | 300 | 600 | 900 | 1200 | 1500 | 2000 | 2500 | 2999 |
|---|---|---|---|---|---|---|---|---|
| clones alive | 98 | 50 | 38 | 27 | 19 | 13 | 11 | 9 |
| times the step | 29,400 | 30,000 | 34,200 | 32,400 | 28,500 | 26,000 | 27,500 | 26,991 |

The sheet obeys both laws of neutral competition, and it had no rule for either: they come from cells that leave by chance and neighbours that grow into the room.

The third result is the one that was not built in. In the sheet, as in the mouse, sisters do opposite things more often than chance: about 55% of pairs had one sister divide and the other leave, where independent sisters would give at most half. Nothing in a cell's network reads its sister, and the networks do not talk to each other. The only thing two sisters share is room: they are born side by side, and when one leaves, the other is the cell most likely to grow into the space and divide. That is the coordination Piedrafita and colleagues proposed, between a cell leaving and a neighbour dividing, and here it is the only coupling there is.

How far below a quarter the mouse's *r* is, and so whether the sheet's sisters are as coordinated as the mouse's, this chapter does not say: it would take the paper's parameter table, which is not used here. Until it is, the sheet's clones should be read as drifting at the mouse's division rate with sister coordination that is real but of unknown strength, and a time read off the sheet in days as an order of magnitude.

## Choices, and why

**Calibrated on the division rate.** It is the one number measured directly, by watching a signal fade, rather than inferred by fitting a model to clone sizes. Calibrating on clone sizes would have assumed the sheet's coordination equals the mouse's.

**In the configuration, not the program.** A tumour cell in the cord divides every twenty steps or so, a basal cell in the sheet every hundred. One step is not one duration for every geometry, so the duration belongs to the configuration that was calibrated, and a configuration that has not been calibrated states none.

**Pairs only once both sisters have acted.** A pair still waiting at the end of a run is left out. That leans towards sisters who acted quickly, and the ten seeds agree to within 0.016 at 1500 steps; at 3000 steps under the default seed *r* is 0.217, so the leaning does not move the answer much.

## Where in the code

`crates/engine/src/sisters.rs`: `SisterFates`, its `r`, and the bookkeeping of lives the loop keeps. `crates/engine/src/simulation.rs`: `Simulation::sister_fates`; where a division begins two lives and a cell leaving ends one. `crates/io/src/config.rs`: `hours_per_step` and `SimConfig::days`. `src/main.rs`: the summary lines, the division interval, the sisters and the scaling of clone sizes. `crates/io/tests/epithelium.rs`: that sisters in the sheet do opposite things more often than chance.

## Still a placeholder

The conversion is the mouse oesophagus's, and only for the sheet: the cord, the spheroid and the stress have no calibrated step. In the human oesophagus cells divide in the first five or six layers rather than in one; nothing here is calibrated against it. Space has no unit but the cell. And the sheet's *r* is compared with the mouse's only by which side of a quarter it falls on.

## References

- Piedrafita G, Kostiou V, Wabik A, et al. [A single-progenitor model as the unifying paradigm of epidermal and esophageal epithelial maintenance in mice](https://doi.org/10.1038/s41467-020-15258-0). *Nature Communications* 11, 1429 (2020). [Free full text](https://europepmc.org/article/PMC/PMC7080751).
