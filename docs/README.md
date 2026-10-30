# The course

What cellweave simulates, explained from the biology up, one notion at a time, in the order the software needed them. Every chapter says what a thing is in the body, what cellweave makes of it, which choices were made and why, where in the code it lives, and what is still a placeholder. The point is to be able to come back in six months and understand the program, not only run it.

Each chapter follows the same shape:

- **In the body**: the biology, defined before it is used.
- **In cellweave**: the model, with a diagram.
- **Choices, and why**: what was decided, and what it costs.
- **Where in the code**: file and symbol, so the claim can be checked.
- **Still a placeholder**: what is a number picked to make something show, not a measurement.
- **References**: the papers the chapter cites, each linked, with a link to its full text where it is free to read.

Chapters build on each other and are meant to be read in order the first time.

| | Chapter | The question it answers |
|---|---|---|
| 1 | [A cell is a set of pixels](01-a-cell-is-a-set-of-pixels.md) | How can a cell have a shape, push on its neighbours, and move, without anyone telling it where to go? |
| 2 | [What happens inside a cell](02-what-happens-inside-a-cell.md) | How does a cell know there is food outside, and decide to grow? |
| 3 | [Oxygen](03-oxygen.md) | Where does oxygen come from, how far does it get, and why is it already at equilibrium whenever a cell looks? |
| 4 | [What a cell does with what it reads](04-what-a-cell-does-with-what-it-reads.md) | Grow, wait, or die: how one number becomes a fate, and how a cell splits in two. |
| 5 | [Two tissues](05-two-tissues.md) | A spheroid that starves from the inside, and a cord that lives around a vessel and never stops turning over. |
| 6 | [A population that evolves](06-a-population-that-evolves.md) | What happens when a daughter is not quite a copy, and nobody says which change is good? |
| 7 | [A tissue that holds its size](07-a-tissue-that-holds-its-size.md) | How does normal tissue stay the same size while renewing itself, and what happens to it as mutations pile up? |
| 8 | [A lethal stress](08-a-lethal-stress.md) | How does a cell die on purpose, and what happens to a tissue when something from outside makes it? |
| 9 | [Rewinding a history](09-rewinding-a-history.md) | Which events made a tumour, and which only came along? |
| 10 | [Many histories, paired](10-many-histories-paired.md) | How often does a hypothesis do what it does, and how are two of them compared fairly? |
| 11 | [Why NOTCH wins](11-why-notch-wins.md) | Which of chapter 7's assumptions let the wrong clones win, and what it takes for NOTCH clones to cover the layer as they do in the oesophagus? |
| 12 | [Units](12-units.md) | How long is a step, and does the sheet obey the laws a real renewing epithelium does? |
| 13 | [Purged by its neighbours](13-purged-by-its-neighbours.md) | What happens to an early tumour that arises among the mutant clones of an ageing epithelium? |
| 14 | [Vessels that grow](14-vessels-that-grow.md) | What happens to a tumour cord when the cells starving at its edge can call new vessels to them? |
| 15 | [A call heard by its neighbours](15-a-call-heard-by-its-neighbours.md) | When a cell calls for vessels, how much of what comes goes to the cells that called? |
| 16 | [Vessels that regress](16-vessels-that-regress.md) | When sprouts that carry no blood are taken back, does a tumour still call in more vessel than it needs? |
| 17 | [Cells pushed aside](17-cells-pushed-aside.md) | How many cells does a growing vessel crush, and do the results of chapter 14 hold when it pushes them aside instead? |
| 18 | [A capillary in the plane](18-a-capillary-in-the-plane.md) | What does a tumour do against a vessel lying across the tissue rather than one seen end on, and can it make blood flow? |
| 19 | [Blood that flows](19-blood-that-flows.md) | When blood goes only where pressure drives it, how much of a tumour's new vessel carries any, and what does the tumour make of it? |
| 20 | [What blood delivers](20-what-blood-delivers.md) | When a vessel delivers oxygen in proportion to the blood through it, what does a tumour make of vessels that carry little? |
| 21 | [Pruned by the flow](21-pruned-by-the-flow.md) | When vessels carrying too little blood are pruned, as real networks prune them, does the network shrink, and does the tumour lose by it? |
| 22 | [Killers that come for the tumour](22-killers-that-come-for-the-tumour.md) | When immune cells arrive through the vessels and kill the tumour cells they touch, do they clear the tumour, hold it, or neither? |
| 23 | [Escape](23-escape.md) | When a tumour's cells can mutate a shield against immune cells, do the immune cells select the shielded ones? |
| 24 | [Measured against the mouse](24-measured-against-the-mouse.md) | Set against the clones counted in the mouse oesophagus, how far apart are the sheet's, and what does it take to bring them together? |
| 25 | [Measured against spheroids](25-measured-against-spheroids.md) | Set against hundreds of measured spheroids, where does the model put the core and the arrest, what does that make of a cord, and how long is a step? |
| 26 | [What the model does not say](26-what-the-model-does-not-say.md) | The placeholders, the missing units, and why emergence is the whole point. |
| | [Glossary](glossary.md) | Every word, in one place. |

## Where the code is

Every chapter ends with where in the code what it describes lives, file and symbol, so that each claim can be checked.

The README at the root of the repository is the front page; this course is the whole account, and its last chapter says what the model does not say.
