# 25. Measured against spheroids

Chapter 5's spheroid was never set against a real one, and its oxygen uptake was chosen to make a core show on a small grid. This chapter sets the model's spheroid against hundreds of measured ones, fixes its uptake by one measured radius, and asks what the rest of the structure, a cord around a vessel, and the model's clock say about it.

## In the body

**Three regions, three radii.** Browning and colleagues (eLife, 2021) grew spheroids from two lines of human melanoma cells, WM983b and WM793b, seeded from about 2,500, 5,000 or 10,000 cells, for up to 24 days. Their cells carry **FUCCI**, the fluorescent ubiquitination cell cycle indicator, two fluorescent proteins that mark a cell in the first gap of its cycle and one past copying its DNA, so that a section through a spheroid shows where cells cycle and where they have stopped; a dye marks the dead. On every section they measured three radii: the spheroid's, the region inside which most cells are arrested, and the necrotic core. They published them for 994 spheroids, under a CC-BY licence.

**Structure follows size, not time.** Their central finding: "comparing spheroid structure as a function of overall size produces results that are relatively insensitive to variability in spheroid size". Spheroids seeded from different numbers of cells pass through the same structures at the same radii: cycling throughout while small; then an arrested centre; then a necrotic core inside the arrested region; and eventually a limiting size, 340 to 356 µm.

**Oxygen sets the core, not all of the arrest.** Murphy and colleagues (PLoS Computational Biology, 2023) grew the same kind of spheroids with their oxygen changed during growth, and found oxygen diffusing in and being consumed enough to account for where the core forms, but not for the arrested region, which a waste product diffusing out accounts for better. Browning and colleagues' model already ascribes arrest to such an inhibitor.

**One radius fixes the oxygen.** With oxygen consumed at one rate by every living cell and held at the surface, a core appears once the centre falls to the level at which cells die, at a radius Murphy and colleagues write $R_c = \sqrt{6 k p_\infty / (\alpha \Omega)}$, from the diffusivity, the level outside and the consumption. Everything about oxygen enters through that one length. In a disc rather than a sphere the 6 becomes a 4, so a disc whose oxygen behaves as the sphere's forms its core at $\sqrt{2/3}$ of the sphere's radius.

**Cords.** Tannock (chapter 14) measured the cords of living tissue around single vessels in a mouse mammary tumour: most had a radius, from the vessel's centre to the necrosis, of 60 to 120 µm.

## In cellweave

**The spheroids, fetched.** `datasets/browning-2021.toml` is the recipe for the authors' published tables (chapter 24 says why the repository keeps recipes and not data), and `scripts/fetch-dataset.py` puts them under `data/`. `scripts/compare-spheroids.py` reads a run's snapshots and works out the same three radii, each as the radius of a disc of the same area: all the tissue; the cells arrested, a hypoxia response at or above the 0.45 at which the simulation stops a cell (chapter 4), or dead; and the necrotic cells. From those, for the runs and the spheroids alike, it gives three numbers:

- $R_c$, where a core appears: the cut between two observed radii that best separates the spheroids with a core from those without;
- $R_q$, where arrest appears, fitted by least squares as Browning and colleagues write the structure before a core forms, $\phi = \sqrt{\max(0, 1 - R_q^2 / R^2)}$, with $\phi$ the arrested radius over the spheroid's;
- the rim of cycling cells, the spheroid's radius less the arrested one, past 1.34 $R_c$, which is where the measured spheroids pass 320 µm.

Runs and spheroids are compared through ratios to their own $R_c$, so that a disc and a sphere are set side by side fairly.

**What the spheroids say.**

| | core from $R_c$ | arrest from $R_q$ | $R_q / R_c$ | cycling rim past 1.34 $R_c$ | rim over $R_c$ |
|---|---|---|---|---|---|
| WM983b, 456 spheroids | 237 µm (97% separated) | 207 µm | 0.87 | 37 µm | 0.16 |
| WM793b, 538 spheroids | 262 µm (96%) | 219 µm | 0.83 | 44 µm | 0.17 |

**The model's thresholds, before anything was run.** At steady state the shipped network arrests a cell below 0.318 of full oxygen, where HIF-1α (hypoxia-inducible factor 1 alpha) reaches 0.45, and kills it below about 0.03, where BCL2 (B-cell lymphoma 2) falls to 0.25. With uptake alone shaping the field, $(R_q/R_c)^2 = (1 - 0.318)/(1 - 0.03)$ in a disc as in a sphere: 0.84.

### The predictions, and what happened

The uptake was set so that the disc's core appears at $\sqrt{2/3}$ of 237 µm, 194 µm, taking a pixel for about a micrometre (chapter 1): chapter 5's formula put death near 50 pixels at 0.0003, and $R_c$ goes as one over the square root of the uptake, which gave 0.00002. Predicted before the run, comparing at the same $R / R_c$: (a) the model's $R_q / R_c$ 0.84, give or take 0.05, above the spheroids'; (b) its cycling rim past 1.34 $R_c$ 62 µm, give or take 15, against their 37; (c) at a fixed radius, its living rim not shrinking over time, where theirs does.

```bash
python3 scripts/fetch-dataset.py datasets/browning-2021.toml
./target/release/cellweave --config examples/spheroid.toml
python3 scripts/compare-spheroids.py --data data/browning-2021/confocal.csv \
  --line 983b output
```

In the first run, under the default seed, the core first appeared between radii of 189 and 204 pixels, a cut of 196 against the 194 aimed at, so the uptake stayed as set; the run was then repeated under seeds 2 and 3. Pooled over the three, read every 25 steps to step 625:

| | core from $R_c$ | arrest from $R_q$ | $R_q / R_c$ | cycling rim past 1.34 $R_c$ | rim over $R_c$ |
|---|---|---|---|---|---|
| the model, 3 runs | 199 px | 168 px | 0.84 | 79 px | 0.40 |

Each run alone gave a core from 196, 201 and 198 pixels, arrest from 168, 168 and 167, a ratio of 0.85, 0.83 and 0.84, and a rim of 80, 79 and 79.

(a) **Held** for the model, 0.84; **failed** against the spheroids, whose ratio first came out at 0.77 and 0.70 from a fit that let spheroids smaller than the arrest radius pull it down, and at 0.87 and 0.83 once fitted as above, before the model was scored. The model puts arrest where the spheroids have it. (b) **Failed** on the number, **held** in direction: the cycling rim is 79 pixels, more than twice the spheroids' relative to their core. (c) **Held**, in a weaker form than stated: the model never stops growing, so no radius stays fixed, but in the first run its living rim stays at 150 to 156 pixels while it grows from 278 to 338, where the spheroids' shrinks from about 140 to 90 µm at a near constant radius.

So oxygen, consumed at one rate and read through the network's two thresholds, gives the core and the onset of arrest at the right sizes: one number fitted, two matched. Late, the spheroids' arrested region reaches much further out than oxygen alone puts it, which is what Murphy and colleagues concluded from changing the oxygen itself.

**A spheroid that never stops.** By step 625, with 802, 827 and 803 dead cells in the three runs, 0, 0 and 3 had been resorbed to nothing. A dead cell's target shrinks by two pixels a step (chapter 4), but it loses a pixel only when a neighbour takes it, and a core of dead cells, each shrinking, has almost no neighbour that will. The core hardly loses volume, and the spheroid grows until the box stops it, where the measured ones stop at 340 to 356 µm. A real core loses mass as its debris breaks down, and Browning and colleagues' model has a rate for it.

### A cord, with the spheroids' oxygen

Predicted before the run: with the uptake fitted to the spheroids, nothing fitted to cords, the cord's radius, the vessel's six pixels and the viable depth, between Tannock's 60 and 120 µm in at least 8 runs of 10.

```bash
./target/release/cellweave --config examples/cord.toml --seeds 1-10 --mcs 1500 --output runs/cord
```

**Held**, in all 10: 68.1 to 70.1 µm at step 1499, 69 in the median, with 50 to 53 living cells. The cord of chapter 5, at its uptake of 0.00003, has about 59. A melanoma in a dish and a mammary tumour in a mouse are not the same cells, so this is a check of order, not a fit: the spheroids' oxygen puts a cord where cords are measured.

### The clock, and what it says

**Space.** Seeded from 2,500 cells, a WM983b spheroid has a radius of 127 µm after a day: 18.7 µm a cell as the diameter of a sphere of the same volume, the space between cells included; 20.4 and 21.9 µm from 5,000 and 10,000 cells. A cell of the model is about 16 pixels across, so a pixel stands for between 1 and about 1.3 µm. The ratios above do not depend on which.

**Time.** Once past a few cells, the model's spheroid doubles about every 50 steps: 32 cells at step 100, 64 at 150, 128 at 200, 256 at 250. Melanoma cells double in 16 to 47 hours, the range Carr and colleagues (Journal of the Royal Society Interface, 2021) gather from two sources, for lines other than these. A step of the tumour is then 0.3 to 0.9 hours.

That reading makes two placeholders measurable:

| what | the model, in those units | measured |
|---|---|---|
| a sprout's tip (chapter 14) | at most a pixel a step: 1.1 to 3.1 µm an hour | 10 to 21 µm an hour in a rabbit's cornea, 5 to 27 in a dish |
| an immune cell's kills (chapter 22) | at most 5 in 300 steps, 4 to 12 days: under 1.25 a day | 1 to 20 a day, the range Weigelin and colleagues cite |
| contact before a kill (chapter 22) | a chance of 0.05 a step: 6 to 19 hours of contact | death after contacts adding up to about an hour or more (Weigelin and colleagues) |

Sprouts are several times too slow, which a model whose tip may not move more than a pixel a step cannot fix without gaps, and they are left as they are. Immune cells kill several times too slowly and take far too long to do it, which two numbers of a configuration can change.

### Immune cells at measured rates

Chapter 22's `cord-immune.toml` is now set to kill with a chance of 0.3 a step of contact, one to three hours of contact at this clock, and at most 20 kills, enough for an immune cell in steady contact to kill within the 1 to 20 a day Weigelin and colleagues cite. When immune cells arrive, how long they live and how strongly they climb the chemokine are not measured, and stay as they were. Predicted before the campaign, at step 1000 against chapter 22's rates: (a) fewer living tumour cells in at least 8 pairs of 10; (b) no living tumour cell left in at least 3; (c) more tumour cells killed in at least 9.

```bash
# cord-immune.toml now kills at the measured rates; chapter 22's, in a copy beside it
sed -e 's/^kill = .*/kill = 0.05/' -e 's/^kills = .*/kills = 5/' examples/cord-immune.toml > examples/cord-immune-22.toml
./target/release/cellweave --config examples/cord-immune-22.toml \
  --versus examples/cord-immune.toml --seeds 1-10 --mcs 1000 --output runs/immune-measured
```

| seed | living, chapter 22's rates | measured rates | killed, chapter 22's | measured |
|---|---|---|---|---|
| 1 | 38 | 0 | 21 | 50 |
| 2 | 37 | 1 | 18 | 43 |
| 3 | 24 | 0 | 61 | 60 |
| 4 | 38 | 0 | 39 | 47 |
| 5 | 33 | 0 | 55 | 53 |
| 6 | 39 | 0 | 22 | 48 |
| 7 | 38 | 0 | 20 | 50 |
| 8 | 18 | 2 | 69 | 47 |
| 9 | 36 | 0 | 41 | 45 |
| 10 | 6 | 0 | 75 | 48 |

(a) **Held**, in all 10. (b) **Held**: 8 cords of 10 had no living tumour cell left, and the other two had one and two; at chapter 22's rates, none. (c) **Failed**: more killed in 6, not 9. A cord cleared early leaves fewer cells to kill, while the slow killers' long histories, the four of chapter 22 in which they held the cord down, killed up to 75.

Chapter 22 found immune cells holding a cord without clearing it. At killing rates several times faster, as measured, they clear most cords within about a thousand steps, two to five weeks at this clock. Whether they would in a real tumour is what the model leaves out: a single kind of immune cell, recognising every tumour cell alike, arriving at a rate nobody measured, against tumour cells that cannot hide (chapter 23 gave them a shield, at chapter 22's rates).

## Choices, and why

**Fitted on the core, not the rim.** One number, the radius at which a core appears, fixes the uptake. Fitting the rim as well would have taken a second free number, and the rim is where oxygen alone is expected to fail.

**The same oxygen as the sphere, not the same radius.** A disc fitted to form its core at the sphere's own 237 µm would have its oxygen reach further than the sphere's does. Fitting the disc at $\sqrt{2/3}$ of it keeps the diffusion length, which is what the cord, where the geometry changes again, can then be checked against.

**Structure against size.** As Browning and colleagues recommend: the model's spheroid starts from one cell and theirs from thousands, and size, not time, is what they share.

**The calibrated values are the defaults.** What a public measurement fixes, the shipped configurations and the program use: the uptake of 0.00002 in every tumour configuration, against a capillary as well as a vessel and in a bath, and as the program's default; the measured killing in `examples/cord-immune.toml` and `examples/cord-evolving-immune.toml`. The chapters whose runs used the values before say so at their top, with the values they used.

## Where in the code

`datasets/browning-2021.toml`, `scripts/fetch-dataset.py`, `scripts/compare-spheroids.py`. `examples/spheroid.toml`, `examples/cord.toml`, `examples/cord-immune.toml`; the program's default uptake in `crates/io/src/config.rs`. The thresholds the comparison reads: `QUIESCENCE_THRESHOLD` and `DEATH_THRESHOLD` in `crates/engine/src/simulation.rs`; the network's in `examples/networks/ras-erk.toml`.

## Still a placeholder

The time unit, from a range of doubling times measured on other melanoma lines. The pixel, between 1 and about 1.3 µm. Three runs of the spheroid, 16 of their snapshots past 1.34 $R_c$. A core that never loses volume. Arrest by oxygen alone, where the spheroids' arrested region grows further than oxygen puts it. Immune cells' arrival, lifetime and pull, unmeasured.

## References

- Browning AP, Sharp JA, Murphy RJ, et al. [Quantitative analysis of tumour spheroid structure](https://doi.org/10.7554/eLife.73020). *eLife* 10, e73020 (2021). [Free full text](https://europepmc.org/article/PMC/PMC8741212).
- Murphy RJ, Gunasingh G, Haass NK, Simpson MJ. [Growth and adaptation mechanisms of tumour spheroids with time-dependent oxygen availability](https://doi.org/10.1371/journal.pcbi.1010833). *PLoS Computational Biology* 19, e1010833 (2023). [Free full text](https://europepmc.org/article/PMC/PMC9876349).
- Carr MJ, Simpson MJ, Drovandi C. [Estimating parameters of a stochastic cell invasion model with fluorescent cell cycle labelling using approximate Bayesian computation](https://doi.org/10.1098/rsif.2021.0362). *Journal of the Royal Society Interface* 18, 20210362 (2021). [Free full text](https://europepmc.org/article/PMC/PMC8455172).
- Weigelin B, den Boer AT, Wagena E, et al. [Cytotoxic T cells are able to efficiently eliminate cancer cells by additive cytotoxicity](https://doi.org/10.1038/s41467-021-25282-3). *Nature Communications* 12, 5217 (2021). [Free full text](https://europepmc.org/article/PMC/PMC8410835).
- Tannock IF. [The relation between cell proliferation and the vascular system in a transplanted mouse mammary tumour](https://doi.org/10.1038/bjc.1968.34). *British Journal of Cancer* 22, 258-273 (1968). [Free full text](https://europepmc.org/article/PMC/PMC2008239).
