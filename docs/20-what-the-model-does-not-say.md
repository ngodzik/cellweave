# 20. What the model does not say

The placeholders, the missing units, and why emergence is the whole point.

## The one rule

Behaviour has to follow from the model, never from a rule that names the behaviour. A cell grows because its ERK is high, not because a flag says grow. A core dies because oxygen does not reach it, not because the code drew a circle. A cord keeps its thickness because births balance deaths, not because a counter holds it there.

This is not a stylistic preference. It is what makes a result mean something. If the necrotic core were placed by a rule, seeing one would say nothing about the hypothesis inside each cell. Because it is not, seeing one at a particular size, on two grids, near where a formula predicts it, says that the network, the diffusion, and the lattice agree with each other and with the closed form. That is a weak claim about cancer and a strong claim about the model, and the second is the one a simulator has to earn first.

The thresholds of chapter 4 are the closest thing to rules: 0.25 on survival, 0.45 on the hypoxia response, 1.4 on volume. They are the points where a continuous signal becomes a discrete fate, and no model avoids them. What they do not do is say where in space a fate lands. That is the test of whether a threshold is a rule in disguise: does it decide a location, or only a level?

## Few units

A pixel is roughly a micrometre, judged by a cell of radius 8 against a real tumour cell of 15 to 20 micrometres across. A lattice step is roughly an hour, judged by a cell dividing in about twenty of them against a real cycle of a day. Roughly, in both cases, and nothing in the code fixes either for the tumour. Until they are fixed, `D = 0.2` cannot be compared to the measured diffusion coefficient of oxygen in tissue, an uptake of `0.0003` cannot be compared to a measured consumption rate, and a cord fifty pixels deep cannot be compared to a cord 150 micrometres deep. Every number that depends on a unit is a placeholder.

**Calibration** is the work of fixing those units and then adjusting the constants until the model reproduces a measured result. For the sheet it has started: chapter 12 fixes a step at about 0.6 hours against the division rate of the mouse oesophagus, and checks two laws of clone sizes and the coordination of sisters. For the tumour it has not: the depth of a real cord, the size at which a real spheroid of a given cell line develops a core and the doubling time of that line are the results to reproduce first, because their numbers are on record.

## The list

What is a placeholder today, gathered from the chapters:

| where | what | chosen so that |
|---|---|---|
| lattice | contact energies 16, 2, 11, 6; λ = 50; T = 10 | cells stay compact and fluctuate visibly |
| network | every rate constant | ERK saturates near 0.55; the network settles in 15 units; survival collapses near 5% oxygen |
| network | growth coefficient 1.5 | a stimulated cell can reach the division threshold at all |
| field | D = 0.2; uptake 0.0003 and 0.00003 | each geometry shows its structure on a 200 pixel grid |
| field | stroma consumes at the tumour's rate | the medium beyond a cord does not fill up from the vessel |
| fates | 0.25, 0.45, 1.4; division checked every 20 steps; resorption 2 pixels per step | the three fates are reachable, in the right order, at the right depths |
| clocks | 8 units of network time per lattice step | the network answers in 2 steps and a cell divides in 20 |
| evolution | mutation rate 1 in 20; factor between 1/2 and 2; every term equally likely | selection acts within a few thousand steps on changes large enough to alter a fate |
| lethal stress | intensity 1; uptake 0.0001; p53 and caspase rates; caspase half way at p53 of 0.4; fragmentation 8 pixels per step | continuous exposure empties the cord against its vessel and leaves a refuge, a periodic one lets it grow back |
| epithelium | NOTCH near 0.01 per step, raised by crowding half way at 60% of reference size; shedding 6 pixels per step; cells 10 pixels across; mutation rate 1 in 100 | the layer turns over in a few hundred steps and holds its size on a 200 pixel grid |
| competition | a ceiling on the size a cell asks for at 1.5 times its reference, what a normal cell asks for; three mutations in ten rewiring | a growth signal buys no size, and a NOTCH term can be lost outright |
| units | 0.6 hours per step for the sheet | a living cell divides about every 97 steps, a mouse oesophageal basal cell every 2.4 days |
| lesion | NOTCH's baseline at 0.7 or 0.3 of its rate; 15 pixels around the centre; at step 4000 | fitter than a normal cell, less fit than a NOTCH-loss clone, about seven cells, in a sheet already covered by mutant clones |
| a silent patch | every living cell within 30 pixels of a point to the right of the vessel, at step 120; read at step 300 | a patch placed before the first sprout, large enough to matter: with the marked network it holds more than half the tumour at step 300 |
| vessels | a wall sprouts at VEGF 0.5; VEGF's reach 30 pixels; tips a pixel per step, 16 pixels apart; new vessel 5 pixels across, perfused at once or once joined, regressing after 100 steps without blood when told to; cells in the way crushed, or pushed along straight lines when told to | a cord whose cells starve calls vessels to them within a few hundred steps, in a 300 pixel box |
| a capillary | seven pixels wide, a quarter of the way down a 300 pixel box; uptake 0.0001; a wall sprouts at VEGF 0.3 | a layer as deep as the cord's sleeve, and walls that sprout at all |
| flow | Poiseuille's law in two dimensions, conductance as the fourth power of the radius over the length; pressure 1 at the inlet, 0 at the outlet; a vessel carries blood above a millionth of the capillary's own flow | blood goes where a difference of pressure sends it, and nowhere else |

## What is not there

- **Only oxygen, crowding and VEGF are coupled.** Growth factor is the same hardcoded input for every cell. VEGF is secreted into a field only where vessels can grow (chapter 14), and what it calls for is a vessel: no cell reads it.
- **Two dimensions.** A tissue here is a monolayer. In three dimensions a sphere starves at a different radius than a disc, and a cell has a dozen neighbours rather than six; the mechanisms are the same, the numbers are not.
- **Tumour and normal tissue meet only in the sheet.** A lesion can be placed in the sheet (chapter 13), as a patch of the basal layer with a network of its own. The tumours of chapters 5 to 10, grown from one cell around a vessel or in a bath, still live in boxes of their own, with no normal tissue around them.
- **Only the field is parallel.** The Monte Carlo step, which moves the pixels, runs on one core. It draws its proposals from one random sequence, and splitting it would change which run a seed names.
- **No picture on disk.** Snapshots record every cell's state and the oxygen it reads, not the grid. The viewer draws the grid live, but saves nothing; the two pictures in chapter 5 came from a debugging dump that is not in the code.

## What a result from this model is

A hypothesis worth testing. Never a result about cancer. The model can say that a given network, under given assumptions about oxygen and mechanics, produces a given tissue, and it can say it reproducibly, on any machine, byte for byte. It cannot say that the network is what is inside a real cell, or that the assumptions hold in a real tumour. Those are for experiments, and the value of the model is in saying which experiment is worth running.

What it can do that an experiment cannot is the counterfactual: the same population, the same seed, one constant changed, so the difference between two runs is attributable to that constant and to nothing else. That is the reason for the determinism check on every push, and it is why every placeholder above is a named constant with a comment rather than a number in the middle of an expression.
