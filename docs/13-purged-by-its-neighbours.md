# 13. Purged by its neighbours

An early tumour arising in normal tissue that is already full of mutant clones, and whether those clones get rid of it.

**Calibrated since.** Chapter 24 set the sheet's NOTCH on the clones counted in the mouse oesophagus: a baseline of 0.001 and crowding at 0.3, half way at 75% of the reference size, where this chapter's runs used 0.005 and 0.03, half way at 60%, and its clones grew about six times too fast per division. A step of the sheet is now about 0.56 hours rather than 0.6. The mechanisms are the same; the shipped examples now run at the calibrated values, and rerunning them gives other numbers than those quoted here.

## In the body

By middle age the lining of the oesophagus is a patchwork of mutant clones (chapters 7 and 11), many of them carrying mutations in genes known from cancer. Yet oesophageal cancer is rare. Colom and colleagues (Nature, 2021) asked what the mutant clones do to a tumour that starts among them.

They followed **early tumours**, microscopic patches of abnormal cells that had only just formed, in mice whose oesophagus had been given many mutations. Most of those tumours were soon gone. Their cells were not dying, not dividing more slowly, and not under attack by the immune system: they were being outcompeted. Making the clones around them more competitive made more tumours disappear; weakening competition made fewer disappear. And the tumours that did survive carried signs of selection. The conclusion: whether an early tumour survives depends on how fit it is compared with the mutant clones around it. The mutant clones of normal tissue, long read as the first steps towards cancer, also get rid of tumours.

What an early oesophageal tumour's cells do differently was measured by lineage tracing (Frede and colleagues, 2016, as reviewed by Alcolea, 2017): a moderate bias towards making dividing cells rather than differentiating ones, and no change in how often they divide.

## In cellweave

**A lesion.** A configuration may place one:

```toml
[lesion]
network = "networks/lesion.toml"   # what its cells carry, next to the configuration
at = 4000                          # the step it appears at
radius = 15.0                      # pixels around the centre of the box
```

At the start of step `at`, every living cell whose centre lies within `radius` pixels of the centre of the box starts carrying the lesion's network, fresh, and all of them join one new clone, marked `"lesion": true` in `clones.json`. Placing it draws nothing, so the run until then is the run without it. A lesion may also say where it is centred, `centre = [x, y]` in pixels, and chapter 15 places one in a tumour.

What makes it a tumour is only its network, which is all the question needs: whether its neighbours outcompete it. `examples/networks/lesion.toml` is the no-pressing epithelium of chapter 11 with NOTCH's baseline slowed to 0.7 of its rate: its cells leave the layer on their own about a third less often than normal cells, and divide no faster, the moderate bias of an early tumour. That makes it fitter than a normal cell and less fit than a cell that has lost NOTCH's baseline altogether, the clones that cover much of the ageing sheet of chapter 11.

**Following it.** `Simulation::lesion_cells` counts the living cells descended from the lesion, in its clone or any clone that arose in it, and `measurements.csv` has it as the `lesion` column, step by step. `Simulation::in_lesion` says whether a given cell is one of them. In the viewer, coloured by clone, the lesion's cells are a pale gold no clone's hue reaches, and the panel follows how many are left.

```mermaid
flowchart LR
    A["the ageing sheet<br/>NOTCH-loss clones spreading"] --> L["step 4000: the cells<br/>within 15 pixels of the centre<br/>carry the lesion's network"]
    L --> N{"fitter than<br/>its neighbours?"}
    N -->|"no"| G["its cells leave the layer<br/>faster than they are made:<br/>lost, no cell killed"]
    N -->|"yes"| S["it holds on, or spreads"]
```

### One history

`examples/epithelium-lesion.toml` is the ageing sheet of chapter 11 in which neither assumption of chapter 7 holds (no pressing, a NOTCH term can be lost, one division in a hundred mutated). The lesion appears at step 4000, about 100 days in at chapter 12's 0.6 hours per step, over six cells, and is followed for 2000 steps more, about 50 days. Mutation goes on, in the lesion too.

Under the default seed:

| step | 4000 | 4200 | 4400 | 4600 | 4800 | 5000 | 6000 |
|---|---|---|---|---|---|---|---|
| lesion cells, ageing sheet | 6 | 7 | 2 | 1 | 1 | 0 | 0 |
| lesion cells, the same sheet with no mutation | 6 | 7 | 11 | 10 | 9 | 7 | 3 |

In the ageing sheet the lesion was gone at step 4866, about 22 days after it appeared. In a sheet that never mutated, the same lesion was still there 50 days on. As in the mice, nothing killed it: no cell died anywhere in the sheet in this history, and every lesion cell that went, left the layer by differentiating, a little faster than the lesion made new ones.

One history is an anecdote, and the two rows are not even the same history: without mutation the generator is drawn differently from the first division on. The next section shows this one was luck.

### Twenty histories, paired

For a fair comparison the clean sheet has to be the same history as the ageing one. A configuration may set `suppress_before`: every mutation founded before that step is drawn exactly as before and not made, the way chapter 9 suppresses one clone. `examples/epithelium-lesion-unmutated.toml` is `epithelium-lesion.toml` with `suppress_before = 4000`. Under one seed the two draw the same numbers: the same divisions, the same departures, the same lesion at the same step, everything but the mutant clones around it. Mutation goes on in both from step 4000. A campaign (chapter 10) now follows the lesion too: the cells left at the end, and the step it was lost at.

```bash
cellweave --config examples/epithelium-lesion.toml --seeds 1-20 --versus examples/epithelium-lesion-unmutated.toml
```

The predictions, written before any code: in the ageing sheet the lesion is lost by step 6000 in at least 12 seeds of 20, in the clean one in at most 6; paired, it ends smaller in the ageing sheet in at least 15; and where it is lost in the ageing sheet, it is lost within 1000 steps (about 25 days) in at least half the seeds, "most were rapidly lost".

**The main ones failed.** The lesion was lost in 16 seeds of 20 in the ageing sheet, and in 16 of 20 in the clean one too. Paired, it ended smaller in the ageing sheet in 3 seeds, larger in 4, and with nothing left in either in 13. Only the timing held: 15 of the 16 losses in the ageing sheet came within 1000 steps, against 10 of 16 in the clean one.

A lesion of six or seven cells that leaves the layer only a third less often than its neighbours is barely fitter than they are, and a clone that small, neutral or nearly so, is lost by chance most of the time whoever its neighbours are: neutral drift, as in chapter 7. So in this model most early tumours of this kind disappear without any competition at all. That is a second explanation for "most were rapidly lost", and Colom and colleagues' other experiments, making the neighbours more or less competitive, are what tells the two apart.

### How fit a lesion has to be

The model can run that experiment directly: keep the two sheets, change the lesion. Written after the first predictions failed and before anything else was run: in the clean sheet the lesion is lost less often as its bias strengthens, at most 8 of 20 at 0.3 of the baseline and at most 4 at none of it; at 0.3, the ageing sheet loses it in at least 4 more seeds than the clean one; and a lesion with no baseline at all, as fit as the fittest clones around it, is not purged by them.

`examples/networks/lesion-strong.toml` is the lesion at 0.3, run by `examples/epithelium-lesion-strong.toml` and `examples/epithelium-lesion-strong-unmutated.toml`; the rows at 0.5 and 0 are `lesion.toml` with NOTCH's baseline rate set to 0.0025 and to 0. Twenty paired seeds each, the sheets as above.

| NOTCH's baseline in the lesion | lost by step 6000, ageing sheet | lost, clean sheet | lesion cells at the end, ageing against clean (medians) |
|---|---|---|---|
| 0.7 of its rate | 16 of 20 | 16 of 20 | 0 against 0 |
| 0.5 | 15 | 2 | 0 against 17 |
| 0.3 | 12 | **0** | 0 against 24.5 |
| none | 0 | 0 | 7 against 85.5 |

All three held. There is a window. A lesion barely fitter than normal cells is lost by drift wherever it lands. A lesion fit enough to hold its own among normal cells, at 0.5 or 0.3, survives in the clean sheet almost every time and is lost in the ageing one most of the time: there, it is the mutant clones that remove it: the two runs of a pair draw the same numbers and differ by those clones and by what they did, nothing else. And a lesion as fit as the fittest clones around it is not removed by them, but it is held back: a median of 7 cells after 50 days among the mutants, 85 in the clean sheet, where it took over more than a quarter of the layer, much as Colom and colleagues (Nature Genetics, 2020) saw equally fit clones stop gaining once they meet.

That is Colom and colleagues' conclusion, reached in the model from the other side: whether an early tumour survives depends on its fitness relative to the clones around it, not on its fitness alone. And it adds what the mice could not show directly: the same lesion, in the same tissue with the same luck, survives or is lost depending only on whether the clones around it carry their mutations.

### Who purged it

Which of those clones? Suppressing a neighbour's mutation from its birth, as chapter 9 does, would not answer it: the mutation was made hundreds of steps before the lesion, and without it the whole sheet would be different by the time the lesion came. So the question is put the way an experiment on mice would put it, by switching the gene back. A scenario may now carry a `Reversion`: at one step, every living cell whose lineage passes through the named clones has those clones' mutations undone, its network rebuilt from the network its lineage started from with every mutation made since but those (`Network` replays a recorded mutation exactly), keeping the protein levels it had. `clones.json` says when, as `reverted_at`. Until that step the run is the real one to the last bit; after it, it differs by those cells and by what they then do.

`cellweave --attribute`, on a configuration with a lesion, runs the real history and notes every clone with a living cell next to a living cell of the lesion while it lived. The candidates are the mutations in those clones' lineages that were acting when the lesion appeared. Each is undone at that step, alone, and then all of them at once, each run on its own core, and for each the table says whether the history stayed identical to the last bit and, if not, what became of the lesion.

The predictions, written before any code, for the strong lesion of the previous section in the first three seeds where it is lost, 1, 2 and 3: undoing all the candidates at once saves it in at least 2 of the 3; in at least 2 of the 3 one single candidate saves it, and every candidate that does carries a NOTCH change that lowers or removes its baseline; and about one candidate per seed leaves the history identical to the last bit, none of them a change to NOTCH, RAS or ERK.

```bash
cellweave --config examples/epithelium-lesion-strong.toml --seed 3 --attribute
```

| seed | the lesion, real history | undoing one candidate | undoing all at once |
|---|---|---|---|
| 1 | lost at step 5139 | 6 of 8 identical to the last bit (three VEGF changes, one BCL2, one HIF-1α, and a link gained by ERK); the two clones that lost NOTCH's baseline, undone one at a time: still lost, earlier, at 4276 and 4446 | 28 cells at step 5999 |
| 2 | lost at step 4863 | two clones, both having lost NOTCH's baseline; undone one at a time: still lost, later, at 4980 and 5735 | 24 cells |
| 3 | lost at step 5746 | clone 436, NOTCH decaying 1.94 times as fast, 52 cells: **saved**, 33 cells at the end; clone 502, NOTCH answering crowding 1.77 times as strongly: saved, 4 cells; 3 of 5 identical (VEGF, BCL2, and a small NOTCH change in a clone of 2 cells) | 29 cells |

The first prediction held: undoing everything that touched the lesion saved it every time. The second failed. Only in seed 3 did one clone decide its fate, a clone keeping NOTCH low; a second rescue there came from a clone answering crowding more strongly, not one with a lower baseline. In seeds 1 and 2 no one clone did it. Each NOTCH-loss neighbour was enough without the others: undo one, and the rest still pushed the lesion out, in seed 1 even sooner than in the real history. The purge was collective, and redundant. The third held on the count, three identical candidates per seed on average, and failed on its clause: a link gained by ERK, and a small NOTCH change in a clone of two cells, were passengers too. Without pressing, a stronger growth signal buys a cell nothing (chapter 11), so an ERK change can do nothing at all.

The passengers are shown, as in chapter 9, rather than inferred: undo them, and every measurement at every step and every pixel at the end come out the same.

### Watching two histories at once

The viewer can run two histories side by side, stepped together batch for batch so that both always show the same step. `--undo` runs, on the right, the same configuration with the named clones undone at the step the lesion appears; `--beside` runs another configuration under the same seed. The panel follows both, the lesion's cells on each side included, and says when the two stop being the same run.

```bash
cd viewer
cargo run --release -- --config ../examples/epithelium-lesion-strong.toml --seed 3 --view clone --undo 436
cargo run --release -- --config ../examples/epithelium-lesion-strong.toml --seed 3 --view clone --beside ../examples/epithelium-lesion-strong-unmutated.toml
```

In the first, the two sheets are the same, pixel for pixel, for 4000 steps. Then the cells of clone 436 turn back into normal cells on the right, and from there the two lesions, pale gold, go their separate ways: lost on the left, at step 5746, still there on the right 50 days on.

## Choices, and why

**A lesion is placed, not evolved.** Waiting for a tumour to arise by mutation in the sheet would make its fitness, its place and its timing different in every history, and the question is about the neighbours. Placing one, with a fitness stated in a file, at a step and a place the configuration fixes, keeps the lesion the same and lets the tissue around it vary.

**Fresh, and one clone.** The lesion's cells start from the lesion network's initial levels, whatever they carried, and share one clone, so that everything descended from it can be followed, and what arose in it told apart from what it started as.

**Placed at the start of a step, drawing nothing.** A lesion that comes after the end of a run leaves the run identical to the one without it, which a test checks.

**A clean sheet by suppression, not by a lower mutation rate.** Running the clean sheet with no mutation would draw different numbers from the first division on, and the two sheets would differ everywhere for reasons that have nothing to do with mutant clones. Drawn and not made, the mutations leave the sheet's luck alone.

**Undone at the lesion's step, not suppressed from birth.** Suppression asks what the tissue would have become without a mutation; undoing it at the lesion's step asks what the clones around the lesion, as they are, do to it. Only the second keeps the lesion's surroundings the same up to the moment it appears, which is what makes a difference afterwards attributable to that clone.

**Candidates from contact.** Colom and colleagues put the competition in the epithelium next to the tumour. A clone that never touched the lesion can still matter, by pushing on a clone that does; it is left out, and undoing all candidates at once saving the lesion every time says the candidates were enough.

**A second set of predictions, said to be one.** The first set failed, and the lesion's bias was then varied. That is the experiment the failure called for, and its predictions were written before it ran; the first lesion is kept as the one the chapter started with.

## Where in the code

`crates/engine/src/simulation.rs`: `Lesion`, the `lesion` field of `Scenario`, and where a step places it; `suppress_before`, read where a division founds a clone. `crates/engine/src/evolution.rs`: the `lesion` field of `CloneRecord`, `Simulation::lesion_cells`, `Simulation::in_lesion`. `crates/engine/src/measure.rs`: the `lesion` column. `crates/io/src/config.rs`: `LesionConfig`, read next to the configuration. `crates/io/tests/lesion.rs`: that a lesion turns the cells at the centre into one clone at its step and not before, that one which never comes changes nothing, that it needs some room inside the box, and that one placed off centre in a tumour takes the cells there.

`crates/engine/src/simulation.rs`: `Reversion`, the `revert` field of `Scenario`, and where a step undoes it. `crates/engine/src/network.rs`: `Network::apply`, which replays a recorded mutation. `crates/engine/src/evolution.rs`: the `reverted_at` field of `CloneRecord`. `src/attribute.rs`: the lesion's attribution, the clones touching it, the candidates and the runs.

`viewer/src/runner.rs`: `Pace`, shared by the runs on screen, `Divergence`, and the dispatch that keeps them in step. `viewer/src/main.rs`: `--undo` and `--beside`.

`src/campaign.rs`: the `lesion` and `lesion_lost_at` columns of `campaign.csv`, and the lesion's line in the summary and the pairs.

`examples/networks/lesion.toml`, `examples/networks/lesion-strong.toml`, `examples/epithelium-lesion.toml`, `examples/epithelium-lesion-unmutated.toml`, `examples/epithelium-lesion-strong.toml`, `examples/epithelium-lesion-strong-unmutated.toml`.

## Still a placeholder

The lesion's bias, four values of it; its size and its place, the centre of the box; the step it appears at, about 100 days into an ageing sheet whose mutation rate is far above a real one. A real early tumour is a three dimensional patch that grows up as well as across; this one is a patch of the basal layer and nothing else. Twenty seeds per row give counts, and where the window begins and ends depends on every placeholder of chapters 7, 11 and 12. Three seeds of attribution, the candidates chosen by contact, and one step to undo them at.

## References

- Alcolea MP. [Oesophageal stem cells and cancer](https://doi.org/10.1007/978-3-319-69194-7_10). A chapter in a Springer book, pages 187-206 (2017). [Free full text](https://europepmc.org/article/PMC/PMC6126191).
- Colom B, Alcolea MP, Piedrafita G, et al. [Spatial competition shapes the dynamic mutational landscape of normal esophageal epithelium](https://doi.org/10.1038/s41588-020-0624-3). *Nature Genetics* 52, 604-614 (2020). [Free full text](https://europepmc.org/article/PMC/PMC7116672).
- Colom B, Herms A, Hall MWJ, et al. [Mutant clones in normal epithelium outcompete and eliminate emerging tumours](https://doi.org/10.1038/s41586-021-03965-7). *Nature* 598, 510-514 (2021). [Free full text](https://europepmc.org/article/PMC/PMC7612642).
- Frede J, Greulich P, Nagy T, et al. [A single dividing cell population with imbalanced fate drives oesophageal tumour growth](https://doi.org/10.1038/ncb3400). *Nature Cell Biology* 18, 967-978 (2016). [Free full text](https://europepmc.org/article/PMC/PMC5870829).
