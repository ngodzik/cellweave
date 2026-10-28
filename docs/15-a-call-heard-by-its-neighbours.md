# 15. A call heard by its neighbours

A cell that calls for vessels calls them to where it is, and its neighbours are there too. This chapter asks how much of what a call brings goes to the cells that made it.

## In the body

Something a cell secretes into the space around it reaches every cell nearby, whether that cell made any or not. Economists call a good that everyone can use and nobody can keep to themselves a **public good**, and the cell that uses one without contributing a **free rider**. Archetti and Pienta (Nature Reviews Cancer, 2019) review what game theory, the mathematics of such situations, says about cancer cells: growth factors that cancer cells produce and that promote the growth of new vessels are among their examples, and "non-cooperative subclones can constantly invade and free-ride on the diffusible factors". What keeps the producers from being overrun is that the good does not go everywhere: "cooperation also evolves more easily if the diffusion range is low", because a producer then keeps more of what it makes, or shares it with producers like itself. How far such factors carry, they add, "remains difficult to estimate".

Axelrod and colleagues (Proceedings of the National Academy of Sciences, 2006) had argued earlier that tumour cells of different lineages can help each other through what they secrete. Salimi Sartakhti and colleagues (PLoS One, 2017) put the call for vessels itself into a spatial public goods game, where what angiogenesis brings goes to neighbours regardless of their contribution.

In a public goods game making the good costs the producer something, and the question is whether producers survive. Here making VEGF costs nothing, which leaves the other half of the question on its own: when a cell calls, who gets the vessel?

## In cellweave

**A patch placed by hand.** The lesion of chapter 13 is a patch of cells that start carrying another network at one step. It may now be placed anywhere, `centre = [x, y]` in pixels, in a tumour as well as a sheet. In a tumour it is a subclone placed by hand, and the `lesion` column of `measurements.csv` follows how many living cells descend from it.

**Louder is not possible.** The first idea was a patch calling twice as loud as the cells around it. It changed nothing: the two runs were the same run to the last bit. Every node of a network is an activated fraction, between zero and one, and a cell starving of oxygen already has its VEGF node at one. In this model a cell cannot call louder than a starving cell does; it can only call less. So the patch calls not at all: `examples/networks/ras-erk-silent.toml` is the tumour network with VEGF never made, and every other node as before.

**Two arms, one history.** `examples/cord-patch-marked.toml` is the cord whose vessels can grow (chapter 14) with a patch placed at step 120, every living cell within 30 pixels of a point to the right of the vessel, whose cells carry the same network as everyone, only marked. `examples/cord-patch-silent.toml` places the same patch with the silent network. Under one seed the two are one run until the patch's silence changes where a vessel goes. Both are read at step 300, and both at three reaches of VEGF, 10, 30 and 100 pixels, the configurations with `reach` changed.

```bash
cellweave --config examples/cord-patch-marked.toml --versus examples/cord-patch-silent.toml --seeds 1-10
```

### The predictions, and what happened

Written before any campaign, after the louder patch had done nothing:

- (a) at a reach of 30 pixels, the silent patch ends with fewer living cells than the marked one in at least 6 pairs of 10: a call brings at least part of what it calls for to the cells that make it;
- (b) at 30 pixels, the rest of the tumour does not lose by it: as many living cells outside the patch or more, with the patch silent, in at least 5;
- (c) the less far a call carries, the more silence costs the patch: the median of the patch's living cells with silence over those without, one added to each, is lowest at 10 pixels and highest, closest to one, at 100.

At a reach of 30 pixels:

| seed | patch, marked | patch, silent | rest of the tumour, marked | rest, silent | sprouts, marked | sprouts, silent |
|---|---|---|---|---|---|---|
| 1 | 2 | 0 | 48 | 45 | 2 | 1 |
| 2 | 31 | 1 | 82 | 76 | 6 | 4 |
| 3 | 16 | 9 | 35 | 39 | 2 | 2 |
| 4 | 41 | 0 | 93 | 39 | 6 | 1 |
| 5 | 35 | 11 | 9 | 38 | 1 | 1 |
| 6 | 55 | 11 | 11 | 27 | 3 | 0 |
| 7 | 52 | 36 | 10 | 1 | 2 | 0 |
| 8 | 41 | 19 | 26 | 19 | 2 | 0 |
| 9 | 27 | 27 | 15 | 15 | 1 | 1 |
| 10 | 37 | 13 | 8 | 24 | 1 | 0 |

(a) **Held**: the silent patch ended smaller in 9 pairs, and the same in the tenth, where the two runs ended alike: there the silence changed nothing the vessels did. In the median it kept about a third of the cells the marked patch did, and was lost entirely in 2 pairs. (b) **Held, just**: the rest of the tumour ended with as many cells or more in 5 pairs, and with fewer in 5. In the pairs where the rest gained (seeds 5, 6 and 10), the silent patch had lost most of its cells; in those where it lost most (seeds 4 and 7), the silent arm had started far fewer sprouts. What a call brings goes first to the cells that make it.

At the three reaches, the patch's living cells with silence over those without, one added to each, and the same for the rest of the tumour (medians over the ten pairs):

| reach of VEGF | runs identical to the last bit | patch, silent over marked | rest, silent over marked | sprouts, marked against silent (medians) |
|---|---|---|---|---|
| 10 pixels | 8 of 10 | 1.00 | 1.00 | 2 against 1.5 |
| 30 pixels | 1 of 10 | 0.35 | 0.97 | 2 against 1 |
| 100 pixels | 10 of 10 | 1.00 | 1.00 | 0 against 0 |

(c) **Failed**, and not by a little: the reach did not trade a private call for a public one. It decided whether a call was heard at all. At 10 pixels the call of cells starving far from the vessel never reached its wall, so silencing them changed nothing, and 8 pairs of 10 were the same run to the last bit. At 100 pixels the call was spread so thin, over tissue that secretes and tissue that does not, that it stayed below the threshold at every wall: no vessel sprouted in any run, with the patch silent or not. Only at 30 pixels did vessels grow and follow who called, and there what a call brought went mostly to the cells that made it.

So in this model a cell that does not call is not a free rider. Vessels come to where the call is strongest, and the cells that called are the ones nearest to it. Whether that holds where VEGF binds to the tissue around it, which shortens its reach and steepens its gradient (chapter 14), is a question for a model that has the binding, and the cost of calling, which this one has neither of.

## Choices, and why

**A patch, not a mutation.** In the evolving cord a mutation that lowers VEGF arises in one cell, at a time and a place chance chooses, and is as likely to die out at once as to found anything. Placing a patch fixes its size, its place and its moment, and leaves only its silence to vary, as the lesion of chapter 13 does for fitness.

**A marked patch to compare with, not none.** The cells of the marked patch carry a fresh copy of the same network, as the silent ones do, so the two arms differ in VEGF and nothing else.

**Silence costs nothing and saves nothing.** A real cell spends something on what it secretes. Here a silent cell divides, arrests and dies exactly as a calling one would in the same place, so whatever the two arms differ by is what the call brought.

## Where in the code

`crates/engine/src/simulation.rs`: the `centre` field of `Lesion`, and where a step places it. `crates/io/src/config.rs`: the `centre` field of `LesionConfig`. `crates/io/tests/lesion.rs`: that a lesion placed off centre in a tumour takes the cells there, and that one outside the box is refused.

`examples/networks/ras-erk-silent.toml`, `examples/cord-patch-marked.toml`, `examples/cord-patch-silent.toml`.

## Still a placeholder

The patch's size, place and step; ten seeds at each reach; three reaches, when how far VEGF carries in a tumour is not known. The cap at one on every node, which is what made a louder call impossible. No cost to calling, so nothing here says whether callers would survive among free riders over many generations.

## References

- Archetti M, Pienta KJ. [Cooperation among cancer cells: applying game theory to cancer](https://doi.org/10.1038/s41568-018-0083-7). *Nature Reviews Cancer* 19, 110-117 (2019). [Free full text](https://europepmc.org/article/PMC/PMC8557269).
- Axelrod R, Axelrod DE, Pienta KJ. [Evolution of cooperation among tumor cells](https://doi.org/10.1073/pnas.0606053103). *Proceedings of the National Academy of Sciences* 103, 13474-13479 (2006). [Free full text](https://europepmc.org/article/PMC/PMC1557388).
- Salimi Sartakhti J, Manshaei MH, Basanta D, et al. [Evolutionary emergence of angiogenesis in avascular tumors using a spatial public goods game](https://doi.org/10.1371/journal.pone.0175063). *PLoS One* 12, e0175063 (2017). [Free full text](https://europepmc.org/article/PMC/PMC5388338).
