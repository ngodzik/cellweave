# 11. Why NOTCH wins

Chapter 7's ageing sheet let the wrong clones win. It named two assumptions that could be to blame; this chapter tests both, together and apart, and only one combination reproduces the oesophagus.

**Calibrated since.** Chapter 24 set the sheet's NOTCH on the clones counted in the mouse oesophagus: a baseline of 0.001 and crowding at 0.3, half way at 75% of the reference size, where this chapter's runs used 0.005 and 0.03, half way at 60%, and its clones grew about six times too fast per division. A step of the sheet is now about 0.56 hours rather than 0.6. The mechanisms are the same; the shipped examples now run at the calibrated values, and rerunning them gives other numbers than those quoted here.

## In the body

Sequencing the oesophagus of people without cancer (Martincorena and colleagues, Science, 2018) found it, by middle age, covered by mutant clones. Two genes stood out. **NOTCH1**, the gene for a receptor that pushes a basal cell to differentiate, was mutated in 12 to 80% of the cells of middle-aged and elderly donors. **TP53**, the gene for p53 (chapter 8), was mutated in 2 to 37%. And NOTCH1 mutations were several times more common in the normal tissue than in oesophageal cancers.

What a NOTCH mutation does to a cell was seen in mice. Alcolea and colleagues (Nature Cell Biology, 2014) switched Notch signalling off in scattered single basal cells and followed their descendants. A normal basal cell divides into two cells that will divide again, two that will differentiate, or one of each, by chance, and on average one of each, so the layer keeps its size (Piedrafita and colleagues, Nature Communications, 2020). In the mutant cells, the divisions that make two differentiating daughters were gone. A mutant clone could no longer be lost by differentiation: it became, in the authors' words, functionally immortal. The mutant cells also pushed their normal neighbours to differentiate, and the mutant clones went on to replace the whole epithelium.

Switching Notch off is a **loss of function**: a change that stops a gene's product from working, rather than turning it up or down a little. It is how Notch is most often hit in the cancers of such epithelia, by inactivating mutations, as Alcolea and colleagues note.

Chapter 7's sheet disagreed with all this. In it, clones that turned the growth pathway up (RAS, the switch that relays a growth signal, and ERK, the relay at the end of it, chapter 2) swept the layer again and again, and clones turning NOTCH down never got past a sixth of it. The chapter put the disagreement down to one of two assumptions: either basal cells do not compete by pressing on each other, or a loss of NOTCH1 does far more than halve a baseline.

## In cellweave

Each assumption can now be switched off on its own.

**Pressing.** In chapter 7 a cell tries to be as large as its growth node says (chapter 2), so a clone with more ERK tries to be larger, squeezes its neighbours, and makes them read crowding and leave. A network's mechanics may now set `max_target`, the largest size a cell asks for, as a multiple of its reference size. `examples/networks/epithelium-no-pressing.toml` is the epithelium network with that ceiling at 1.5, which is what a normal cell asks for. A stronger growth signal then buys no size and no pressure; a weaker one still lowers the size a cell asks for, so a cell short of growth signal still cannot divide. A ceiling below the division ratio of 1.4 would stop every cell dividing, and is refused.

**Losing a term.** Chapter 6 gave mutations the power to rewire: with `rewiring` set, a share of mutations lose a term outright or gain a link. NOTCH has three terms: a slow pull towards full on, its baseline, which makes a cell leave now and then even with room around it; a decay; and the rise with crowding. A scaled rate can at most halve or double one of them. A lost baseline is a loss of function: a cell that only leaves when its neighbours squeeze it.

**What sequencing measures.** `Simulation::coverage` gives, for every node, the share of the living cells whose lineage carries a change to it, whatever the change does. That is what a sequencing study reports gene by gene. A run with mutation prints it at the end; a campaign writes it as the `coverage` column of `campaign.csv` and summarises it per arm.

Four configurations make the two by two:

| | mutations only scale a rate | a mutation can lose a term |
|---|---|---|
| **cells press by growing** | A: `epithelium-mutating.toml` (chapter 7) | B: `epithelium-rewiring.toml` |
| **no cell presses harder than a normal one** | C: `epithelium-no-pressing.toml` | D: `epithelium-no-pressing-rewiring.toml` |

```bash
cellweave --config examples/epithelium-mutating.toml --seeds 1-10 --versus examples/epithelium-rewiring.toml
cellweave --config examples/epithelium-no-pressing.toml --seeds 1-10 --versus examples/epithelium-no-pressing-rewiring.toml
```

Each is the ageing sheet of chapter 7: 6000 steps, one division in a hundred mutated, seeds 1 to 10. Each pair takes about two minutes on twenty threads.

### The predictions, and what happened

Written before any code, comparing at the end the share of the layer carrying a NOTCH change with the share carrying a change to the growth pathway (the larger of the RAS and ERK shares):

- A: growth ahead of NOTCH in at least 8 seeds of 10. **Held, in 10.**
- B: growth still ahead in at least 6, because a cell that stops leaving on its own is still pressed out by a neighbour that grows. **Held, in 10.**
- C: NOTCH ahead in at least 8, its share between 12 and 80%. **Held: 8 seeds, a median of 19%.**
- D: NOTCH ahead in at least 9, a median share of at least a half, and a lost NOTCH term in the leading line in at least 5. **Held: 10 seeds, a median of 63%, a lost term in 6.**

| | mutations only scale a rate | a mutation can lose a term |
|---|---|---|
| **cells press by growing** | A: every living cell, in every seed, carries a change to RAS and to ERK; NOTCH 13% (median) | B: the same, every cell carrying growth changes; NOTCH 36%, its baseline lost in 4 leading lines, each of them already turned up on growth |
| **no cell presses harder than a normal one** | C: no sweep; NOTCH 19%, growth 2% | D: NOTCH 63% (19 to 99%), growth 2% (none above 14%); the leading line lost NOTCH's baseline in 6 seeds |

Pressing decides everything when it is allowed. In A and B the growth pathway took the whole layer in every seed, and a NOTCH loss, when it came, came on top, in a line that had already won by size.

Without pressing and without losses, in C, nothing sweeps. NOTCH's 19% is barely more than the 14% of BCL2 (B-cell lymphoma 2, the survival protein of chapter 2), whose changes can do nothing in a sheet at full oxygen and spread only by drift. Halving a chance of about one in a hundred per step is too small an advantage to show in 6000 steps.

Only D looks like the oesophagus. NOTCH changes hold more of the layer than growth changes in every seed, between 19 and 99% of it, within the 12 to 80% sequenced in people in 6 seeds and above it in 4; growth changes stay rare. And the loss that wins is always the same one: in each of the 6 leading lines that lost a NOTCH term, it is the baseline, never the decay and never the response to crowding. That is the mutant of Alcolea and colleagues: a clone that is no longer lost by leaving on its own, and still leaves when squeezed.

So in this model the oesophagus's pattern needs both: basal cells that do not compete by size, and NOTCH mutations that abolish what they hit rather than dim it. Either alone is not enough. That is a claim about the model, and a hypothesis about the tissue: that a growth signal does not buy a basal cell room.

## Choices, and why

**A ceiling, not a new growth law.** One number in the network file is the smallest change that removes pressing, and the rest of the mechanics stays as chapters 2 and 4 describe it. Set at what a normal cell asks for, so that it takes nothing from a normal cell and gives nothing to one with more growth signal.

**Coverage counts any change.** A sequencing study counts mutations in a gene whether or not they do anything, and a change that does nothing spreads only with what it travels with. Comparing NOTCH with a node whose changes cannot matter, as BCL2 in C, is what separates selection from drift.

**Two campaigns, not one with four arms.** A campaign compares two arms. Within each pair the two arms draw the same numbers until the first mutation, where rewiring draws one more and the histories part, so the comparison between arms is of outcomes, not of pairs.

## Where in the code

`crates/core/src/network.rs`: `max_target` on `MechanicsSpec`. `crates/engine/src/network.rs`: the ceiling applied where the growth node becomes a target size; `Network::node_names`. `crates/engine/src/simulation.rs`: the refusal of a ceiling below `DIVISION_RATIO`. `crates/engine/src/evolution.rs`: `Simulation::coverage`. `src/campaign.rs`: the `coverage` column and the per-arm summary. `src/main.rs`: the shares printed after a run.

`examples/networks/epithelium-no-pressing.toml`, `examples/epithelium-rewiring.toml`, `examples/epithelium-no-pressing.toml`, `examples/epithelium-no-pressing-rewiring.toml`.

## Still a placeholder

Everything chapter 7 lists, the ceiling at exactly the normal size, three mutations in ten rewiring, and ten seeds of each. Two things the papers saw are not measured yet: whether a NOTCH-loss clone makes its normal neighbours leave faster, as Alcolea and colleagues saw, and whether two such clones slow when they meet, as Colom and colleagues did (Nature Genetics, 2020). TP53 has no counterpart here: nothing damages a cell in the sheet, so a cell that ignores damage has no advantage. And the steps are still steps; what 6000 of them are in the life of a mouse or a person is the next chapter's question.

## References

- Alcolea MP, Greulich P, Wabik A, et al. [Differentiation imbalance in single oesophageal progenitor cells causes clonal immortalization and field change](https://doi.org/10.1038/ncb2963). *Nature Cell Biology* 16, 615-622 (2014). [Free full text](https://europepmc.org/article/PMC/PMC4085550).
- Colom B, Alcolea MP, Piedrafita G, et al. [Spatial competition shapes the dynamic mutational landscape of normal esophageal epithelium](https://doi.org/10.1038/s41588-020-0624-3). *Nature Genetics* 52, 604-614 (2020). [Free full text](https://europepmc.org/article/PMC/PMC7116672).
- Martincorena I, Fowler JC, Wabik A, et al. [Somatic mutant clones colonize the human esophagus with age](https://doi.org/10.1126/science.aau3879). *Science* 362, 911-917 (2018). [Free full text](https://europepmc.org/article/PMC/PMC6298579).
- Piedrafita G, Kostiou V, Wabik A, et al. [A single-progenitor model as the unifying paradigm of epidermal and esophageal epithelial maintenance in mice](https://doi.org/10.1038/s41467-020-15258-0). *Nature Communications* 11, 1429 (2020). [Free full text](https://europepmc.org/article/PMC/PMC7080751).
