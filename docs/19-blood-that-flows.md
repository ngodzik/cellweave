# 19. Blood that flows

Until now a vessel either carried blood or did not, by a rule: from its first pixel, or once its tip had joined another vessel. This chapter lets pressure decide, as it does in the body, and asks how much of a tumour's new vessel ends up carrying any.

**Calibrated since.** Chapter 25 fixed the tumour's oxygen uptake on measured melanoma spheroids at 0.00002, where this chapter's runs used 0.0001. The mechanisms are the same; the shipped examples now run at the calibrated value, and rerunning them gives other numbers than those quoted here.

## In the body

**Pressure drives blood.** Blood runs from where its pressure is high to where it is low, through every path open to it, and each path takes a share that depends on how easily blood passes through it. For a narrow tube that is **Poiseuille's law**: the flow is the difference of pressure across the tube times the fourth power of its radius, over its length and the blood's viscosity. A tube twice as wide carries sixteen times as much; one twice as long, half.

**Vessels that carry no flow.** In a network, a vessel open at both ends may still carry nothing, if its ends stand at the same pressure, and a blind end never carries anything. Owen and colleagues (Journal of Mathematical Biology, 2009) made blood flow through a model network of vessels growing in tissue, and found pruning by low flow "highly sensitive to the pressure drop across a vascular network". Pries and colleagues (American Journal of Physiology, 1998) compared simulated networks with every segment of real ones in the rat mesentery, the sheet of tissue that holds the gut, estimating the flow through each segment with a model of flow through the network.

## In cellweave

**A network of tubes.** With a capillary across the box (chapter 18), a configuration may say:

```toml
[angiogenesis]
perfusion = "flow"   # the others are "at_once" and "on_joining"
```

Then the vessels are a network: junctions where a sprout starts on a vessel, where a tip has joined one, and where a blind tip ends, joined by segments along the capillary and along each sprout. Each segment's conductance is the fourth power of its radius over its length, three pixels for the capillary and two for a sprout, the constant every tube shares dropped. Blood enters at the capillary's left end, held at pressure 1, and leaves at its right, held at 0. At every other junction, what flows in flows out, which makes the pressures a linear system; `crates/engine/src/flow.rs` solves it exactly, by elimination, and gives each segment its flow. A sprout carries blood, and so oxygen, where the flow through it is more than a millionth of what the capillary alone would carry.

The flow is solved again whenever the network changes in a way that can change it: when a tip joins a vessel, and when a sprout regresses (chapter 16). A sprout that the flow leaves is dry from then, and with `regress_after` it regresses that many steps later.

**What this changes.** A blind sprout carries nothing, as with `on_joining`. But a loop now carries blood only if there is a difference of pressure across it. A sprout that starts on a blind sprout and comes back to it closes a loop with no pressure across it, and carries nothing. A loop between two nearby points of the capillary carries little, and one that bypasses a long stretch of it carries a lot.

**Watching it.** In the viewer a vessel is as bright as the blood through it: full red at the capillary's own flow, a third darker for each tenfold less, dark where none flows. Coloured by `flow`, living cells are a plain grey, dead ones near black, and the network over them runs from faint violet to white. `Simulation::flow_through` says what flows through a vessel pixel and `Simulation::flowing` whether blood flows as pressure drives it at all.

```bash
cd viewer
cargo run --release -- --config ../examples/capillary-delivery.toml --view flow
```

`examples/capillary-joining.toml` is `capillary-angiogenic.toml` with blood only through loops; `examples/capillary-flow.toml`, with blood where pressure drives it.

### The predictions, and what happened

Written before the campaign, after one run of the default seed in which 67% of the vessel carried blood by flow against 90% once joined: at step 400, (a) a smaller share of the vessel carrying blood by flow than once joined, in at least 9 pairs of 10; (b) fewer living cells by flow in at least 7; (c) more necrotic cells by flow in at least 7.

```bash
cellweave --config examples/capillary-joining.toml --versus examples/capillary-flow.toml --seeds 1-10 --mcs 400
```

| seed | living, once joined | by flow | necrotic, once joined | by flow | vessel carrying blood, once joined | by flow | sprouts, once joined | by flow |
|---|---|---|---|---|---|---|---|---|
| 1 | 164 | 70 | 1 | 12 | 82% | 70% | 49 | 29 |
| 2 | 153 | 98 | 0 | 9 | 92% | 75% | 46 | 38 |
| 3 | 139 | 83 | 2 | 10 | 87% | 93% | 31 | 26 |
| 4 | 153 | 106 | 0 | 13 | 83% | 64% | 40 | 38 |
| 5 | 148 | 93 | 1 | 11 | 77% | 68% | 43 | 33 |
| 6 | 119 | 79 | 2 | 19 | 85% | 62% | 40 | 40 |
| 7 | 135 | 89 | 3 | 11 | 76% | 72% | 42 | 38 |
| 8 | 128 | 101 | 1 | 11 | 83% | 94% | 42 | 33 |
| 9 | 143 | 99 | 2 | 20 | 76% | 79% | 40 | 42 |
| 10 | 173 | 112 | 0 | 15 | 82% | 82% | 47 | 32 |

(a) **Failed**: a smaller share by flow in 7 pairs, not 9; 73% against 83% in the median. (b) **Held**, in all 10: 95.5 living cells against 145.5. (c) **Held**, in all 10: 11.5 necrotic cells against 1.

The share of vessel carrying blood moved little, and the tumour still lost a third of its living cells. (These are the numbers of the run after a review of the code found a join onto a sprout placed one disc early in the network; the first run, before the fix, had come to the same three outcomes.) What moved was how much vessel there was: by flow, the tumour started fewer sprouts, 35.5 against 42, laid less vessel, 5838 pixels against 7154, and so had about a third less vessel carrying blood. A loop that closes feeds no one unless pressure drives blood round it, and the tumour fed by flow was fed less and starved more. That it also started fewer sprouts was not predicted, and why was not tested.

## Choices, and why

**Two dimensions, three-dimensional tubes.** The network lies in the plane, but each segment conducts as a round tube of its radius. A two-dimensional channel would conduct as the cube of its width; the tube is what a capillary is.

**Pressure held at the two ends of the capillary.** What drives blood through the window is the difference of pressure between where the capillary enters and where it leaves. Only ratios of flow matter to who is fed, so the pressures are 1 and 0, and the flows are read as shares of what the capillary alone would carry.

**Exact, and only when needed.** A growing network has at most a few hundred junctions, so the linear system is solved exactly, by elimination with partial pivoting, rather than iterated. A tip growing on a blind sprout changes no flow, so the system is solved again only when a tip joins a vessel or a sprout regresses.

**Carrying blood, or not.** Above a millionth of the capillary's own flow, a vessel is held at full oxygen like any other; below, it is held at nothing. How much oxygen a vessel gives off should depend on how much blood runs through it, which is the next thing to add.

## Where in the code

`crates/engine/src/flow.rs`: `Segment`, `conductance`, `pressures`, `flows`. `crates/engine/src/vessels.rs`: `Capillary`, `Join`, where a tip records the vessel it joined, `Tree::network`, which builds the segments and junctions, and `Tree::solve_flow`, which shares the flow out over each sprout's discs and says which sprouts carry blood and since when they have been dry. `crates/core/src/types.rs`: `Perfusion::Flow`. `crates/engine/src/simulation.rs`: where a flow is refused without a capillary.

`viewer/src/paint.rs`: the `Flow` view, `flow_strength` and `vessel_colour`. `viewer/src/glow.rs`: a vessel's light by its flow.

`examples/capillary-joining.toml`, `examples/capillary-flow.toml`.

## Still a placeholder

Blood is a fluid of one viscosity, with no red cells to share out unevenly at a branch. Every sprout has one radius, which never adapts to the flow it carries. A vessel carrying blood holds full oxygen whatever its flow. The pressures at the two ends are fixed whatever the network downstream does.

## References

- Owen MR, Alarcón T, Maini PK, Byrne HM. [Angiogenesis and vascular remodelling in normal and cancerous tissues](https://doi.org/10.1007/s00285-008-0213-z). *Journal of Mathematical Biology* 58, 689-721 (2009).
- Pries AR, Secomb TW, Gaehtgens P. [Structural adaptation and stability of microvascular networks: theory and simulations](https://doi.org/10.1152/ajpheart.1998.275.2.H349). *American Journal of Physiology, Heart and Circulatory Physiology* 275, H349-H360 (1998).
