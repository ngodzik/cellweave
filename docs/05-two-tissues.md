# 5. Two tissues

A spheroid that starves from the inside, and a cord that lives around a vessel and never stops turning over.

Everything in chapters 1 to 4 is a mechanism. This chapter is what the mechanisms do together, in two geometries that the configuration chooses between. Both are classic experimental systems, and both produce a structure nobody drew.

## The spheroid

### In the body, or rather in the dish

Grow tumour cells in a drop of liquid so they cannot attach to anything, and they clump into a ball: a **spheroid**. It is the standard way to study a tumour in three dimensions without an animal. Small spheroids are alive throughout. Past a few hundred micrometres across, the centre dies: a **necrotic core**, ringed by quiescent cells, ringed by a proliferating rim. The layering is set by how far oxygen and nutrients diffuse in before being used up.

### In cellweave

One cell in the middle of a bath. It grows, divides, and its daughters do the same. Every cell reads the oxygen where it sits, and as the ball grows, the cells in the middle read less. Past a certain size they stop cycling, and past a larger one they die.

```
#################xx##################################
################xxxxxxxxx###########xx###############
##########xxx##xxxxxxxxx#######xxxxx#################
########x#xxxxxx#xxxxxxxx##xx#xxxxxxxx###############
#########xx##xxxxxxxxxxxxxxxxxxxxxxxxx###############
###########xxxxxxxxxxxxxxxxxxxxxxxxxxxxx#############
#########xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx#x#######
######xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx#########
########xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx#x#x#######
#########xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx#######
#########xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx#xxxxx#
##xxx####xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx##
#xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx####
###xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx#####
##xxxxxx##xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx#x######
###x######xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx#######
###########xxxxxxxxxxxxxxxxxxxxxxxxxxx##xxxxxx#######
########xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx##xxxxxxx######
#######xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx#####xxx#x#####
######xxxxxxx###xxxxxxxxxxxxxxxxxxxxx#########xx#####
######xxxxxx#x#######xxxxxxxxxxxxxx##################
######################xxxxxxxxxxxxx##################
#############################xxx#####################
##############################xx#x###################

#  living     x  necrotic          200 x 200 grid, step 200, seed 42
```

The core is not placed. It appears because oxygen does not get that far, and because the network turns a lack of oxygen into a collapse of survival. And it appears at the same size on two grids one and a half times apart, which is what the bath is for:

| | 200 x 200 | 300 x 300 |
|---|---|---|
| first necrotic cell | step 150 | step 150 |
| cells at that moment | 59 | 58 |
| equivalent radius | 61 px | 60 px |

The bathed disc formula of chapter 3 predicts the centre reaching the death level near a radius of 50; the first death comes a little later, which is what a network that needs two steps to answer, on a tissue that is not a perfect disc, should give.

### Its limit

The rim always has oxygen, so it always proliferates, so the ball always grows outward. On the default grid it touches the edge of the box at about step 216, and the run stops there, since past that point the box would be deciding what the tissue sees. Ten or so generations. And nothing evolves in a spheroid: cells reproduce only at the rim, where selection is nil, and die at the core, where they do not reproduce. A tolerant cell in the core dies later and passes nothing on.

## The cord

### In the body

Cut a lung cancer thin and look at it under a microscope, as was first done in the 1950s, and the living tissue is arranged in sleeves around each blood vessel, with necrosis filling the space between the sleeves. Each sleeve is a **tumour cord**: alive out to the distance oxygen diffuses from the vessel, dead beyond. This is the founding observation of tumour hypoxia. It is why "necrotic beyond about 150 micrometres from a vessel" is the best known number of tumour hypoxia, and it is a steady state: cells are born near the vessel, pushed outward, and die at the edge, indefinitely.

### In cellweave

An obstacle disc in the middle, held at full oxygen: the vessel. Stroma all around, consuming. One cell seeded against the wall.

```
..xxxxxxxx###################################xxxxxxxxxxxxxx.
xxxxxxxxxx##################################xxxxxxxxxx#xxxxx
.xxxxx#xxx####################################xx#x#xxx#xxxxx
xxxx#####################################################xxx
xxxx#####################################################xxx
xxxx####################################################xxxx
xxxxxx#x###################################################x
xxxxx######################################################x
xxxxx####################################################x#x
xxxxx######################################################x
xxxx#######################################################x
xxxx######################################################xx
xxxx####################################################xxxx
xx######################################################xxxx
x#############################O##########################xxx
x###########################OOOOO#########################xx
xx#########################OOOOOOO#########################x
x###########################OOOOO#########################xx
xx############################O#############################
xxxx########################################################
xxxx########################################################
xx#########################################################x
x###########################################################
xx#########################################################x
xx#########################################################x
xxxx#xx####################################################x
xxxxxx#####################################################x
xxxxxxxxxx##################################################
xxxxxxxxxx################################################xx
xxxxxxxxxxx##############################################xxx
xxxxxxxxxxx############################################xxxxx
..xxxxxxxxx#########################################xxxxxxxx

O  vessel    #  living    x  necrotic, being resorbed    .  medium
200 x 200 grid, step 800, seed 42
```

The geometry is the spheroid's inverted: living inside, dead outside. And because the dead are resorbed, the tissue does not grow. Cells are born near the wall, are pushed outward, stop cycling, die at the viable boundary, and are cleared, and the ring keeps its thickness.

| | value |
|---|---|
| living cells, steps 350 to 3000 | 41 to 45 |
| born by step 2950 | 1052 |
| gone by step 2950 | 963 |
| turnover | one birth and one death every three steps |
| edge of the box | never touched |

A flux balance at the vessel wall predicted about 43 living cells. The first prediction, made with the spheroid's formula, was off by a factor of two, and that is recorded in the experiment log: a vessel is a small source whose flux thins with distance, not a bath.

### Why it matters

This is the tissue an evolving population needs. Every birth is a chance to vary. A lineage that tolerates less oxygen survives further out; if it still divides there, its daughters inherit that, and the cord thickens. Thousands of generations in a stable environment, on a grid the tissue never fills. Chapter 6 turns mutation on and watches what selection does with it.

## What both taught that the prototype hid

Three things showed up only when the mechanisms were run together, and each looked like biology until it was checked:

- With one unit of network time per lattice step, a cell's fate trailed its surroundings by most of a generation, and the spheroid's core appeared at 112 cells rather than a dozen.
- Continuing a spheroid run after the tissue touched the box left the bath nowhere to be; every cell then read near zero and all 158 died, rim included.
- Starting the cord's field from zero left a slow front that the solver's tolerance never caught. The first cell drifted from 0.28 to 0.32 over 1400 steps, crossed the quiescence threshold, and the whole population appeared at step 1500. Nothing biological happened at step 1500.

## How a tissue is measured

Both tissues are judged on the same few numbers, which every run writes to `measurements.csv`, one row per step, next to its snapshots:

- **Cycling, quiescent, necrotic, differentiating**: how many cells are in each fate of chapter 4, and how many are leaving the layer of an epithelium (chapter 7).
- **Born, gone**: cells created and cells resorbed since the start. The difference between two rows is the turnover over that window.
- **Clones**: how many clones still have a living cell. One, unless mutation is on (chapter 6).
- **Tissue radius**: the radius of a disc with the area the tissue covers, vessel included.
- **Necrotic radius**: the same for the necrotic cells alone. In a spheroid this is the core.
- **Viable depth**: how far living tissue reaches from where its oxygen comes in. From a bath, the tissue radius less the necrotic radius, since the living are a ring outside the core. From a vessel, the radius of vessel and living together less the vessel's own, since the living are a sleeve around it.

Every radius is an **equivalent radius**, the radius of a disc of the same area, because a tissue is never a perfect disc and this is the one radius that does not depend on where its boundary wobbles. It is also why the numbers are rough on a tissue far from round: a ragged rim covers more area than a smooth one of the same depth.

Read off a run, these are the numbers this chapter quotes. The cord's viable depth holds between 52 and 55 pixels from step 500 to step 3000, and over its second half it makes 38.5 cells and resorbs 39.0 per hundred steps. The spheroid, when it stops at the edge, has a necrotic radius of 53 pixels and a viable depth of 38.

The run itself ends with the same in three lines, the turnover taken over the second half of the run so that the growth from a single cell does not count.

## Where in the code

`crates/engine/src/measure.rs`: `Measurements`, and `Simulation::measure`, which reads the lattice and the fates and draws nothing, so measuring a run does not change it. `crates/io/src/output.rs`: `MeasurementLog`, the table.

`crates/core/src/types.rs`: `OxygenSource`, bath or vessel with a radius. `crates/io/src/config.rs`: `oxygen_uptake`, `seed`, and `SimConfig::scenario`, which turns a configuration into what the engine runs. The command line's `--seed` and `--output` run one configuration under another seed, somewhere else, without editing it.

`examples/spheroid.toml` and `examples/cord.toml`: the two configurations, run twice each by CI and compared byte for byte.

`crates/engine/src/simulation.rs`: `declare_sources` is where the geometry decides who consumes and who supplies; `touches_the_edge` is the stop.

The experiment log in loomkeeper holds the two experiments with their hypotheses written before the runs, and the corrections made along the way.

## Still a placeholder

The spheroid's viable rim is about one cell thick where a real one keeps several; the cord's depth of about fifty pixels against a vessel of radius six is a ratio of eight where the body's is nearer thirty. Both follow from uptake values chosen to make the structure visible on a 200 pixel grid. Neither has been compared to a measured spheroid or a measured cord, and doing so is what calibration will mean.
