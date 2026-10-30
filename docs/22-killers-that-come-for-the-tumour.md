# 22. Killers that come for the tumour

Every tumour of this course has grown in tissue with no immune system. This chapter lets killer immune cells arrive through the vessels, find the tumour and kill the cells they touch, and asks whether they clear it, hold it, or neither.

**Calibrated since.** Chapter 25 fixed the tumour's oxygen uptake on measured melanoma spheroids at 0.00002, where this chapter's runs used 0.00003. The mechanisms are the same; the shipped examples now run at the calibrated value, and rerunning them gives other numbers than those quoted here. Chapter 25 also set immune cells' killing on measured rates, a chance of 0.3 a step of contact and at most 20 kills, where this chapter's runs used 0.05 and 5: at those rates immune cells clear 8 cords of 10.

## In the body

**Killers that patrol.** A **killer T cell**, or cytotoxic T lymphocyte, recognises a cell that is infected or abnormal and kills it by contact, punching holes in its membrane and setting off its apoptosis (chapter 8). Killer T cells leave the blood through the walls of small vessels and crawl through tissue, drawn by **chemokines**, small proteins cells secrete that attract other cells up their gradient, which is **chemotaxis**.

**Slower killers than once thought.** Halle and colleagues (Immunity, 2016) watched killer T cells at work in living mice, against cells infected by a virus, and found that "on average, one CTL killed 2-16 virus-infected cells per day", that they "remained migratory" while killing, and that a cell contacted by more than two of them was likelier to die. Weigelin and colleagues (Nature Communications, 2021) found that in solid tumours killer T cells "often fail to kill target cells during 1:1 conjugation", and that a tumour cell died when several of them added sublethal damage, "3 serial hits delivered with decay intervals below 50 min".

**Elimination, equilibrium, escape.** Dunn and colleagues (Annual Review of Immunology, 2004) named what an immune system does to a tumour **immunoediting**, in three phases: "elimination, equilibrium, and escape". It may clear the tumour; it may hold it, alive but contained; or the tumour may come to escape it, the immune system "sculpting the immunogenic phenotype of tumors as they develop" by killing the cells it sees best.

## In cellweave

**Immune cells.** A configuration with a vessel may add:

```toml
[immunity]
from = 200          # the step immune cells start to arrive
arrival = 0.1       # chance per step that one arrives
kill = 0.05         # chance per step that one kills a tumour cell it touches
kills = 5           # tumour cells one can kill before it is exhausted
lifetime = 300      # steps one lives
chemotaxis = 200.0  # how strongly it climbs the chemokine
reach = 20.0        # pixels over which the chemokine falls by e
```

From step `from`, each step, with chance `arrival`, one immune cell arrives: a small cell of the lattice, 13 pixels, a new kind, `CellKind::Immune`, placed on a pixel beside a vessel wall carrying blood, drawn among all such pixels. It carries no network, so it does not grow, divide or call for vessels, and it is not counted among the living. `CpmLattice::seed_cell` places it, taking pixels from whatever cells held them as a growing vessel does.

**It climbs the chemokine.** Every living tumour cell secretes a chemokine into a field of its own, which diffuses and decays like VEGF (chapter 14), falling by e over `reach`. In the Cellular Potts step an immune cell extending one pixel up the chemokine's gradient saves `chemotaxis` times the difference, `CpmLattice::attract_immune`; it sticks to tumour cells, less than they stick to each other, and keeps away from other immune cells. It keeps to its size loosely, a crawling lymphocyte changing shape as it goes.

**It kills what it touches.** Each step, an immune cell touching a living tumour cell kills the first of them in reading order with chance `kill`: the tumour cell starts apoptosis and breaks up as in chapter 8. After `kills` kills, or `lifetime` steps, it dies and breaks up in turn. Arrivals and kills are drawn from the run's generator, and only with immunity, so every run without it is what it was.

**What a run says.** `measurements.csv` gains `immune`, the immune cells in the tissue, and `killed`, the tumour cells they have killed; a run's summary prints both, and a campaign writes them and says in how many runs no living tumour cell was left. In the viewer immune cells are an icy white, and the panel counts them.

`examples/cord-immune.toml` is `cord-wide.toml` (chapter 14) with immune cells from step 200.

### A first version that could not move

The first immune cells kept to their size as firmly as a tumour cell does. For a cell of 13 pixels that held it still, every pixel gained or lost costing more than the chemokine's pull could pay, and the immune cells killed only where they had arrived, beside the vessel: from step 200 the cord fell from 42 living cells to 11 by step 800, with 83 cells killed. The predictions first written were written for that version; it was fixed before any campaign, and the predictions written again.

### The predictions, and what happened

Written again after fixing that, before the campaign, from one run of the default seed in which the motile immune cells killed 30 tumour cells by step 900 and the cord held 46 living cells against about 40 without them: at step 1000, with immune cells, (a) fewer living cells in at least 6 pairs of 10; (b) no living cell left in at most 1; (c) among the living, a larger share quiescent in at least 7, the killers arriving through the vessel and meeting the well-fed cycling cells first.

```bash
cellweave --config examples/cord-wide.toml --versus examples/cord-immune.toml --seeds 1-10 --mcs 1000
```

| seed | living, without | with immune cells | cycling, without | with | quiescent, without | with | necrotic, without | with | immune cells at the end | tumour cells killed |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 40 | 38 | 8 | 8 | 32 | 30 | 43 | 41 | 0 | 21 |
| 2 | 37 | 37 | 8 | 8 | 29 | 29 | 49 | 47 | 1 | 18 |
| 3 | 40 | 24 | 9 | 2 | 31 | 22 | 44 | 8 | 11 | 61 |
| 4 | 38 | 38 | 8 | 5 | 30 | 33 | 44 | 40 | 1 | 39 |
| 5 | 37 | 33 | 7 | 3 | 30 | 30 | 44 | 19 | 8 | 55 |
| 6 | 38 | 39 | 7 | 8 | 31 | 31 | 47 | 40 | 0 | 22 |
| 7 | 37 | 38 | 8 | 8 | 29 | 30 | 45 | 40 | 0 | 20 |
| 8 | 38 | 18 | 8 | 2 | 30 | 16 | 45 | 0 | 19 | 69 |
| 9 | 40 | 36 | 7 | 8 | 33 | 28 | 45 | 32 | 3 | 41 |
| 10 | 39 | 6 | 8 | 0 | 31 | 6 | 45 | 0 | 25 | 75 |

(a) **Held**, just: fewer living cells in 6 pairs, more in 2, as many in 2; 36.5 against 38 in the median. (b) **Held**: no tumour was cleared, though seed 10 was down to 6 living cells. (c) **Failed**: a larger share quiescent in 6, not 7.

The immune cells held the cord in two kinds of history. In four (seeds 3, 5, 8 and 10) they killed 55 to 75 tumour cells and the cord fell to between 6 and 33 living cells; in most of the others they killed about 20, and the cord was the cord without them. In every pair there were fewer necrotic cells with immune cells, 36 against 45 in the median, and none at all in seeds 8 and 10. Two readings fit, and neither was tested: a tumour cell killed breaks up four times as fast as one that starves (chapters 4 and 8), and one killed leaves its oxygen to the cells beyond it. Nor was what decides which kind of history a cord has. Few immune cells were left at the end, 2 in the median: most of those that came had used their kills, lived their time, or been squeezed out by the tumour cells around them.

These are the numbers of the run after a review of the code. Its first run had counted immune cells squeezed out of every pixel as still there, about thirty at the end of each run; it had also let a dying immune cell at the edge of the box count as tissue reaching it, and two immune cells kill one tumour cell in the same step.

## Choices, and why

**Through vessels carrying blood.** Killer T cells reach tissue from the blood, so they arrive where blood flows, at a place drawn among every pixel beside such a wall.

**From a step on.** An immune response takes days to start, and a tumour of one cell beside the vessel immune cells arrive through is killed at once (the first version did exactly that). Starting at step 200 lets a cord form first.

**A chance per step, a limited count, a limited life.** A kill chance of one in twenty per step of contact, with five kills over 300 steps, makes one immune cell kill about two cells a day on the reading of a step as some thirteen minutes (chapter 14), at the low end of Halle and colleagues' 2 to 16. The cooperation and the added hits Weigelin and colleagues saw are not modelled: one hit kills.

**A chemokine from every living tumour cell.** What draws killer T cells into a real tumour is made by tumour cells, immune cells and the tissue around them alike, in response to much. Here every living tumour cell secretes it alike, which is the least that brings immune cells to the tumour at all.

## Where in the code

`crates/core/src/types.rs`: `Immunity`, `CellKind::Immune`. `crates/engine/src/cpm.rs`: `CpmLattice::seed_cell`, `CpmLattice::attract_immune`, the chemotaxis term in a flip, and the contact energies with immune cells. `crates/engine/src/simulation.rs`: `IMMUNE_LAMBDA`, where immune cells arrive, kill and die, the chemokine, `Simulation::immune_cells`, `Simulation::killed`, `Simulation::chemokine`. `crates/engine/src/measure.rs`: the `immune` and `killed` columns. `src/campaign.rs`, `src/main.rs`: where they are reported. `viewer/src/paint.rs`: `IMMUNE`. `viewer/src/panel.rs`: the immune line.

`examples/cord-immune.toml`.

## Still a placeholder

Every number of the `[immunity]` table. Chapter 25 sets two of them against measured killing, at a clock fixed on spheroids: there they are several times too slow, and at the measured rates immune cells clear 8 cords of 10 by step 1000. One kind of immune cell, recognising every tumour cell alike: no antigen, no priming in a lymph node, no exhaustion but a count. One hit kills. The chemokine comes from tumour cells alone.

## References

- Dunn GP, Old LJ, Schreiber RD. [The three Es of cancer immunoediting](https://doi.org/10.1146/annurev.immunol.22.012703.104803). *Annual Review of Immunology* 22, 329-360 (2004).
- Halle S, Keyser KA, Stahl FR, et al. [In vivo killing capacity of cytotoxic T cells is limited and involves dynamic interactions and T cell cooperativity](https://doi.org/10.1016/j.immuni.2016.01.010). *Immunity* 44, 233-245 (2016). [Free full text](https://europepmc.org/article/PMC/PMC4846978).
- Weigelin B, den Boer AT, Wagena E, et al. [Cytotoxic T cells are able to efficiently eliminate cancer cells by additive cytotoxicity](https://doi.org/10.1038/s41467-021-25282-3). *Nature Communications* 12, 5217 (2021). [Free full text](https://europepmc.org/article/PMC/PMC8410835).
