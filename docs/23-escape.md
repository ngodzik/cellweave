# 23. Escape

A tumour whose cells mutate is a population under selection (chapter 6). With killer immune cells in the tissue (chapter 22), a change that shields a cell from them is one selection could favour. This chapter gives the tumour's network such a shield, lets it mutate, and asks whether immune cells select the cells that carry more of it.

## In the body

**A shield on the cell's surface.** **PD-L1**, programmed death-ligand 1, is a protein on a cell's surface that binds PD-1, programmed cell death protein 1, on a killer T cell, and holds the killer back. Many tissues show it, which keeps an immune response from doing more damage than it must; many tumours show it too. Juneja and colleagues (2017) deleted PD-L1 from mouse tumour cells and found that PD-L1 on the tumour cells "potently inhibited CD8+ T cell cytotoxicity", CD8 marking the killer T cells, and that tumour cells keeping it "outcompeted PD-L1-deleted MC38 cells in vivo, demonstrating tumor PD-L1 confers a selective advantage".

**Escape, the third E.** Dunn and colleagues (Annual Review of Immunology, 2004) describe the immune system "sculpting the immunogenic phenotype of tumors as they develop": by killing the tumour cells it recognises best, it leaves the others, and a tumour that started out visible can end up escaping it. That is selection, with the immune system as what selects.

## In cellweave

**A shield in the network.** A network may name an **evasion** role, the node whose level shields the cell: an immune cell touching it kills it with the chance it would otherwise, times one less that level. `examples/networks/ras-erk-pdl1.toml` is the tumour network of `ras-erk.toml` with a PD-L1 node, switched on by the growth factor every cell reads and turned over, which settles near 0.17. A mutation scaling either of its rates by up to two moves it to between about 0.09 and 0.29; and like any node's, its rates mutate only when a daughter mutates, one of every rate of every node equally likely (chapter 6). `Network::evasion` reads it; `Simulation::evasion` is its mean over the living tumour cells, which a campaign writes and compares pair by pair.

`examples/cord-evolving-pdl1.toml` is the cord of chapter 14 in its wide box, with that network and a mutation in one division of twenty; `examples/cord-evolving-immune.toml` adds the immune cells of chapter 22. The two are one run under one seed until the first immune cell arrives at step 200.

```bash
cellweave --config examples/cord-evolving-pdl1.toml --versus examples/cord-evolving-immune.toml --seeds 1-10 --mcs 2000
```

### The predictions, and what happened

Written before anything was run with the PD-L1 network, at step 2000: (a) a higher mean PD-L1 among the living with immune cells than without, in at least 7 pairs of 10; (b) with immune cells, above its starting 0.17, in at least 7; (c) without immune cells, within 0.02 of 0.17 in at least 7, drift alone moving it little. A run stops when its tumour reaches the edge of the box.

| seed | living, without immune cells | with | PD-L1 among the living, without | with | tumour cells killed | reached the edge, with immune cells |
|---|---|---|---|---|---|---|
| 1 | 39 | 42 | 0.167 | 0.166 | 57 | no |
| 2 | 62 | 170 | 0.131 | 0.167 | 26 | at step 703 |
| 3 | 37 | 96 | 0.167 | 0.187 | 48 | at step 1971 |
| 4 | 39 | 226 | 0.118 | 0.167 | 26 | at step 1122 |
| 5 | 47 | 39 | 0.171 | 0.165 | 68 | no |
| 6 | 72 | 167 | 0.182 | 0.167 | 65 | at step 1908 |
| 7 | 40 | 52 | 0.167 | 0.197 | 40 | no |
| 8 | 42 | 41 | 0.168 | 0.167 | 73 | no |
| 9 | 36 | 187 | 0.109 | 0.167 | 40 | at step 1350 |
| 10 | 32 | 9 | 0.167 | 0.167 | 91 | no |

(a) **Failed**: a higher mean PD-L1 with immune cells in 5 pairs only. (b) **Failed**: above 0.17 in only 2, seeds 3 and 7. (c) **Held**: within 0.02 of 0.17 in 7 without immune cells.

The immune cells did not select the shield. A mutation raising PD-L1 to 0.29 cuts a cell's chance of being killed by about a seventh, and the immune cells killed some 26 to 91 tumour cells in 2000 steps; in most runs no mutation of PD-L1 spread among the living either way.

The tumours escaped all the same, and not that way. Without immune cells no tumour reached the edge of the box in 2000 steps. With them, five did, between steps 703 and 1971, with 96 to 226 living cells against a cord of about forty. In two of them every living cell carried a change to HIF-1α (hypoxia-inducible factor 1 alpha, chapter 2), the kind of change by which a line escapes its hypoxic arrest and grows out of the cord (chapter 10); but in seed 9 no mutation at all was shared by the living, so mutation is not the whole of it. Why immune cells should make a tumour outgrow its cord was not tested. One reading is that every kill frees room and oxygen for a division, so that a tumour under attack divides more, and with mutation finds an escape sooner; another, that a cord thinned by killing reorganises around its vessel in a way that grows.

## Choices, and why

**A shield, not a switch.** PD-L1 holds a killer back more or less, so its level scales the chance of a kill rather than deciding it.

**Made by the growth factor.** Real tumour cells raise PD-L1 when signals from the immune response reach them, as well as by mutation. Here it answers only the constant growth factor, so that what moves it is mutation alone, and anything that rises has been selected.

**No cost.** A cell making more PD-L1 pays nothing for it here, so whatever changes its share of the tumour is the immune cells.

## Where in the code

`crates/core/src/network.rs`: the `evasion` field of `Roles`. `crates/engine/src/network.rs`: `Network::evasion`. `crates/engine/src/simulation.rs`: where a kill is shielded, and `Simulation::evasion`. `src/campaign.rs`: the `evasion` column and its tally.

`examples/networks/ras-erk-pdl1.toml`, `examples/cord-evolving-pdl1.toml`, `examples/cord-evolving-immune.toml`.

## Still a placeholder

PD-L1's rates and its baseline; that it shields in proportion to its level; that it costs nothing; that nothing but mutation raises it. Only one way to escape is open, where real tumours also stop showing what killers recognise.

## References

- Dunn GP, Old LJ, Schreiber RD. [The three Es of cancer immunoediting](https://doi.org/10.1146/annurev.immunol.22.012703.104803). *Annual Review of Immunology* 22, 329-360 (2004).
- Juneja VR, McGuire KA, Manguso RT, et al. [PD-L1 on tumor cells is sufficient for immune evasion in immunogenic tumors and inhibits CD8 T cell cytotoxicity](https://doi.org/10.1084/jem.20160801). 214, 895-904 (2017). [Free full text](https://europepmc.org/article/PMC/PMC5379970).
