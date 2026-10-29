# 21. Pruned by the flow

In chapter 16 a sprout regressed only if it carried no blood at all. Real networks prune more than that: they remove the segments that carry little. This chapter prunes by flow, and asks whether a tumour's network shrinks, and whether the tumour loses by it.

## In the body

**Low flow, pruned.** Franco and colleagues (PLoS Biology, 2015), whom chapter 16 quotes, proposed that flow "leads to stabilization of high flow/shear segments and regression of adjacent low flow/shear segments", the shear being the drag blood exerts on the vessel's wall, which grows with the flow. Owen and colleagues (Journal of Mathematical Biology, 2009) built pruning by low shear into a model of tissue whose vessels grow and carry blood, and found it "highly sensitive to the pressure drop across a vascular network". In both, what decides whether a segment stays is how much blood runs through it, not merely whether any does.

## In cellweave

**A floor on the flow.** With blood flowing as pressure drives it (chapter 19) and sprouts that regress (chapter 16), a configuration may add:

```toml
[angiogenesis]
perfusion = "flow"
regress_after = 100
prune_below = 0.02
```

A sprout whose flow stays below `prune_below`, as a share of what the capillary alone carries, counts as dry: it regresses `regress_after` steps after its flow fell below the floor, with every sprout grown from it, as one carrying no blood does. Every sprout of one radius, the shear on its wall goes as its flow, so a floor on the flow is a floor on the shear. Without `prune_below`, only a sprout carrying no blood at all regresses, as before.

`examples/capillary-regressing.toml` is `capillary-delivery.toml` (chapter 20) with `regress_after = 100`; `examples/capillary-pruning.toml` adds `prune_below = 0.02`, a fiftieth of the capillary's own flow, a placeholder.

### The predictions, and what happened

Written before anything was run with it, at step 600: with pruning, (a) more sprouts regressed in at least 8 pairs of 10; (b) less vessel at the end in at least 7; (c) the living cells within a tenth of the unpruned twin's in at least 6, pruning removing vessels that deliver little.

```bash
cellweave --config examples/capillary-regressing.toml --versus examples/capillary-pruning.toml --seeds 1-10 --mcs 600
```

| seed | living, unpruned | pruned | vessel pixels, unpruned | pruned | regressed, unpruned | pruned | sprouts, unpruned | pruned |
|---|---|---|---|---|---|---|---|---|
| 1 | 64 | 71 | 8417 | 9427 | 43 | 52 | 99 | 112 |
| 2 | 88 | 88 | 10094 | 10318 | 31 | 43 | 96 | 107 |
| 3 | 90 | 75 | 9374 | 8051 | 21 | 30 | 76 | 80 |
| 4 | 70 | 80 | 8790 | 9989 | 12 | 23 | 65 | 83 |
| 5 | 73 | 83 | 8623 | 9508 | 41 | 50 | 99 | 112 |
| 6 | 95 | 98 | 9285 | 9177 | 17 | 23 | 78 | 81 |
| 7 | 85 | 88 | 8498 | 8376 | 36 | 37 | 91 | 88 |
| 8 | 89 | 90 | 8723 | 8424 | 30 | 43 | 88 | 99 |
| 9 | 86 | 83 | 8272 | 7622 | 30 | 37 | 84 | 86 |
| 10 | 91 | 80 | 9658 | 8241 | 23 | 36 | 85 | 91 |

(a) **Held**, in all 10: 37 sprouts regressed against 30 in the median. (b) **Failed**: less vessel in only 6 pairs, more in 4; 8800 pixels against 8757 in the median, the same. (c) **Failed**: the living cells stayed within a tenth of the unpruned twin's in only 5 pairs. In the others pruning moved them by up to a sixth, more in some pairs and fewer in others: more living cells in 6 pairs, fewer in 3, as many in 1, and 83 against 87 in the median.

Pruning removed more, and the network did not shrink. Every pruned sprout freed the wall it had started from, which sprouted again while the cells around it still called, as chapter 16 found for blind ends: the pruned tumours started 89.5 sprouts against 86.5. The network was renewed, not reduced. Which vessels it ended with, and so how many cells lived, changed with the history, in no one direction.

These are the numbers of the run after a review of the code found three faults in what this chapter builds on: a join onto a sprout placed one disc early, weak vessels held as sinks of oxygen, and a weak sprout pruned with a strong one grown from it. The first run, before the fixes, had found pruning cost the tumour nothing in all ten pairs; that did not survive them.

## Choices, and why

**The flow, not the shear.** The shear on a vessel's wall goes as its flow over the cube of its radius. With one radius for every sprout, a floor on one is a floor on the other, and the flow is what the model already has.

**The same delay as for dry sprouts.** A sprout carrying too little regresses after the same `regress_after` as one carrying none, counted from when its flow fell below the floor. One delay keeps the comparison with chapter 16 to one difference.

**A sprout at a time.** A sprout is pruned whole, with what grew from it, where a real network removes one segment between two junctions. A sprout carrying blood through its first half and none beyond keeps its blind half.

## Where in the code

`crates/core/src/types.rs`: the `prune_below` field of `Angiogenesis`. `crates/engine/src/vessels.rs`: where solving the flow marks a sprout strong or weak, and since when, and where a weak sprout regresses. `crates/engine/src/simulation.rs`: where `prune_below` is refused without flow and regression.

`examples/capillary-regressing.toml`, `examples/capillary-pruning.toml`.

## Still a placeholder

The floor, a fiftieth of the capillary's flow. Radii that never adapt: in a real network a segment carrying much blood widens and one carrying little narrows, and Pries and colleagues (chapter 19) needed several stimuli besides the shear to keep a simulated network from collapsing to a few paths. Whole sprouts pruned, not segments.

## References

- Franco CA, Jones ML, Bernabeu MO, et al. [Dynamic endothelial cell rearrangements drive developmental vessel regression](https://doi.org/10.1371/journal.pbio.1002125). *PLoS Biology* 13, e1002125 (2015). [Free full text](https://europepmc.org/article/PMC/PMC4401640).
- Owen MR, Alarcón T, Maini PK, Byrne HM. [Angiogenesis and vascular remodelling in normal and cancerous tissues](https://doi.org/10.1007/s00285-008-0213-z). *Journal of Mathematical Biology* 58, 689-721 (2009).
- Pries AR, Secomb TW, Gaehtgens P. [Structural adaptation and stability of microvascular networks: theory and simulations](https://doi.org/10.1152/ajpheart.1998.275.2.H349). *American Journal of Physiology, Heart and Circulatory Physiology* 275, H349-H360 (1998).
