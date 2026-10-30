# 20. What blood delivers

In chapter 19 a vessel carrying any blood at all held full oxygen at its wall, whether blood rushed through it or barely moved. This chapter makes what a vessel delivers follow the blood through it, and asks what a tumour makes of vessels that carry little.

**Calibrated since.** Chapter 25 fixed the tumour's oxygen uptake on measured melanoma spheroids at 0.00002, where this chapter's runs used 0.0001. The mechanisms are the same; the shipped examples now run at the calibrated value, and rerunning them gives other numbers than those quoted here.

## In the body

**Oxygen comes with the blood.** Oxygen travels bound to the red cells of the blood, and a vessel gives the tissue around it only what the blood brings. Where blood barely moves, little arrives, and what arrives is soon given off; where it rushes, the vessel stays rich along its whole length.

**A tumour's blood flow is uneven.** Netti and colleagues (Microvascular Research, 1996) open their study of how pressure and flow are related in tumours with: "Tumor blood flow (TBF) is characterized by spatial and temporal heterogeneities". Some of a tumour's vessels carry plenty, some little, and which ones changes.

## In cellweave

**Delivery that follows flow.** With blood flowing as pressure drives it (chapter 19), a configuration may add:

```toml
[angiogenesis]
perfusion = "flow"
half_delivery_flow = 0.1
```

A vessel then delivers a share of what it would at full flow: nothing without blood, half at a flow of `half_delivery_flow`, as a share of what the capillary alone carries, and all of it at the capillary's own flow, $\frac{q (1 + h)}{q + h}$ for a flow $q$ and half $h$, never more than all. It gives and never takes: where the tissue around a vessel pixel already holds more than it would deliver, the pixel is left free, as a vessel carrying no blood is, rather than held low; held low, a vessel carrying a trickle would draw oxygen out of the tissue as if its blood could carry any amount away. The capillary delivers by the same rule, column by column: where a bypass takes part of its blood, the stretch it bypasses carries less. Each vessel pixel is held at that share of full oxygen, and of the lethal stress (chapter 8) when an exposure is on, since the stress arrives with the blood too. `ScalarField::hold_pixels_by` holds each pixel at a value of its own, and `hold_vessels_by` holds the vessels at what they deliver.

`examples/capillary-delivery.toml` is `capillary-flow.toml` with `half_delivery_flow = 0.1`.

### Exploring first

One run of the default seed, to step 400, at two halves. At a hundredth of the capillary's flow, little changed: 100 living cells against 93 delivering in full, most sprouts carrying far more than a hundredth. At a tenth, the tumour had 71 living cells, started 52 sprouts against 45 and laid 7797 pixels of vessel against 6491. The example uses a tenth.

### The predictions, and what happened

Written before the campaign, after that run: at step 400, with delivery following flow, (a) fewer living cells in at least 7 pairs of 10; (b) more sprouts started in at least 7; (c) more vessel and fewer living cells, both, in at least 6.

```bash
cellweave --config examples/capillary-flow.toml --versus examples/capillary-delivery.toml --seeds 1-10 --mcs 400
```

| seed | living, full delivery | following flow | necrotic, full | following flow | vessel pixels, full | following flow | sprouts, full | following flow |
|---|---|---|---|---|---|---|---|---|
| 1 | 70 | 68 | 12 | 15 | 5514 | 6001 | 29 | 30 |
| 2 | 98 | 75 | 9 | 17 | 6793 | 7954 | 38 | 52 |
| 3 | 83 | 67 | 10 | 21 | 4713 | 4962 | 26 | 32 |
| 4 | 106 | 64 | 13 | 23 | 6538 | 7585 | 38 | 49 |
| 5 | 93 | 83 | 11 | 18 | 5231 | 5922 | 33 | 37 |
| 6 | 79 | 65 | 19 | 16 | 6311 | 6623 | 40 | 45 |
| 7 | 89 | 65 | 11 | 15 | 6222 | 6598 | 38 | 41 |
| 8 | 101 | 78 | 11 | 21 | 4680 | 5395 | 33 | 37 |
| 9 | 99 | 75 | 20 | 19 | 6162 | 6440 | 42 | 42 |
| 10 | 112 | 76 | 15 | 18 | 5488 | 6606 | 32 | 43 |

(a) **Held**, in all 10: 71.5 living cells against 95.5 in the median. (b) **Held**, in 9: 41.5 sprouts against 35.5, as many in the tenth. (c) **Held**, in all 10: every tumour whose vessels delivered by flow had more vessel and fewer living cells than its twin, 6519 pixels against 5838. (These are the numbers of the run after a review of the code; the first run, before it, had come to the same three outcomes.)

More vessel, and less fed. Vessels carrying a little blood deliver a little oxygen, so the cells around them go on starving and calling, more walls sprout, and the new vessels, carrying little themselves, deliver little in turn. Nothing in the model says a tumour should have many vessels and poor ones; the call does what it always did, and the vessels it brings give less.

## Choices, and why

**A share that saturates.** What a vessel can give rises with the blood through it and cannot exceed what full blood gives. The curve chosen is the simplest that does both, with one number, the flow at which it gives half.

**Anchored on the capillary.** The capillary alone, carrying its own flow, delivers in full, so that a run with no sprouts is the same whatever the half: before any vessel grows, nothing changes.

**The stress comes with the blood.** Whatever a vessel brings, it brings in proportion to its blood, so the lethal stress is held at the same share as the oxygen.

## Where in the code

`crates/core/src/types.rs`: the `half_delivery_flow` field of `Angiogenesis`. `crates/engine/src/pde.rs`: `ScalarField::hold_pixels_by`. `crates/engine/src/coupling.rs`: `hold_vessels_by`. `crates/engine/src/vessels.rs`: the flow along the capillary column by column, and `Vessels::flow_through`. `crates/engine/src/simulation.rs`: `delivered`, and where the oxygen and the stress are held by it.

`examples/capillary-delivery.toml`.

## Still a placeholder

The half, a tenth of the capillary's flow. A vessel delivers the same share along its whole length, where real blood gives off its oxygen as it goes and arrives poorer at the far end. Blood has no red cells to share out unevenly where a vessel branches. The flow, once solved, does not fluctuate.

## References

- Netti PA, Roberge S, Boucher Y, et al. [Effect of transvascular fluid exchange on pressure-flow relationship in tumors: a proposed mechanism for tumor blood flow heterogeneity](https://doi.org/10.1006/mvre.1996.0041). *Microvascular Research* 52, 27-46 (1996).
