# 10. Many histories, paired

One history proves little. Twenty, run side by side and compared two by two under the same luck, show which outcomes a hypothesis allows and how often, and here they overturned both predictions made for them.

## In the body

The same kind of tumour, in different people, does different things, and part of the reason is history: which mutations happened to arise, in which cell, at which moment (chapter 9). Faced with the same exposure, one population is wiped out and another comes back. So a question about how a population responds to something is not answered by one history but by many, and the honest answer is a distribution: in how many histories this, in how many that.

## In cellweave

`--seeds 1-20` runs one configuration under every seed of a range, all at once, one seed per core. `--versus other.toml` adds a second configuration, an **arm**, run under the same seeds. Each run reports where it ended: the step, whether it touched the edge, the cells in each fate, the clones alive, the fewest living cells from the start of an exposure on, the share of the living that is resistant, their mean response to the stress, and the latest mutations of the clone leading at the end. `campaign.csv` holds one row per seed and arm, and the run prints a summary: each arm on its own, then the two compared seed by seed.

Seed by seed, because of **common random numbers**. Two arms that differ only by their exposure schedule are, under the same seed, the same history until the step the schedules first differ: the same divisions, the same mutations, the same tissue. Comparing them pair by pair takes that shared history out of the comparison, so what differs between the two runs of a pair is what differs between the arms and not the luck of the draw. It is a known method in the study of simulations, and the exact reproducibility of a run (chapters 1 and 9) is what makes it available here.

**Resistance** is measured on each living cell's network, out of the tissue: given full oxygen and the stress at half its source's level, about what the cells nearest the source read, what is its chance per step of starting apoptosis, compared with the unmutated network's? One is unchanged; below one is less sensitive; below a half counts as resistant. Asked out of the tissue so that a cell sheltered in a refuge is not taken for resistant because the stress does not reach it.

### Continuous or intermittent

`examples/cord-evolving-stress.toml` is the evolving cord of chapter 6, mutating one division in twenty, with the stress response network of chapter 8, exposed from step 800 continuously. `examples/cord-evolving-stress-periodic.toml` differs by one line: the stress is on for 150 steps in every 500. Seeds 1 to 20, 3000 steps each, paired.

The predictions, written before any code: the intermittent arm ends with more living cells in at least 15 pairs of 20; and resistant cells make up a larger share of the living under the continuous arm in at least 12. **Both were wrong.** The intermittent arm had more living cells in 12 pairs (fewer in 7, one tie), and a larger resistant share went to the continuous arm in 6 pairs and to the intermittent in 4, the other 10 tied, almost all at none.

One change was made between the predictions and the campaign, and it is recorded here because it could have been made to fit: the resistance measure was first set at the full level of the source, where the response saturates and no single mutation can halve it. A trial on two seeds counted nothing as resistant, and the measure was moved to half the level before the twenty seeds were run. The predictions were not touched.

What happened instead is more interesting than either prediction. Under continuous stress the outcome is all or nothing:

| continuous, 20 seeds | how | seeds |
|---|---|---|
| **collapsed** to a refuge, 1 to 17 cells, none dividing | the stress empties the cord against its vessel (chapter 8) | 2, 3, 4, 6, 8, 9, 10, 12, 14, 15, 19 |
| **survived by resistance**, the living responding at 27 to 62 percent of the unmutated network | a mutation turning p53 or caspase down in the leading line | 5, 7, 13, 18 |
| **survived without resistance**, the living responding as the unmutated network does | the cord escapes its hypoxic arrest and grows out of the stress's reach | 1, 11, 16, 17, 20 |

The third route was not expected, and was checked rather than assumed. In seeds 11 and 20, at step 800 the cycling cells read enough stress to raise p53 to a median of 0.19 and 0.16; at the end every living cell is cycling and reads almost none, p53 at a median of 0.01 and 0.02. The stress is taken up faster than oxygen (its uptake is more than three times oxygen's), so it reaches less far. A line that has escaped arrest, as the lines of chapter 6 do, keeps dividing out past where the stress reaches, into tissue the oxygen still feeds. Every one of those five leading lines carries a change to HIF-1α, and in seed 20 the leading line is the same in both arms: it was on its way before the exposure began.

Under intermittent stress the cord never collapses: it ends with 26 to 90 living cells in every seed, touching the edge in 7. A line that has turned its response down, to under three quarters of the unmutated network's, leads in 5 seeds of 20 (5, 11, 14, 16, 18), as one does under continuous stress (5, 7, 13, 18, and 9 among the collapsed), but mostly not in the same seeds: the two arms of a pair are one history until step 950, when the intermittent stress first goes off, and part from there.

So, in this model, continuous stress controls the cord far better in about half of the histories and fails in the others, and which histories it fails in is decided by mutations that were there before it started or that arose against it; intermittent stress never eliminates the cord, and does not make resistance rarer. That is a statement about this network, this geometry and these rates, and nothing more; but it is the kind of statement, how often and by which route, that only many paired histories can make.

## Choices, and why

**Pairs, not means.** Twenty runs of each arm compared as two clouds would mix the difference between the arms with the difference between histories. Paired by seed, the histories are the same up to the exposure, and the comparison is of the arms.

**Medians and counts, not averages.** The outcomes are bimodal, collapse or not, and an average of 1 and 88 describes neither. The summary gives medians and ranges, and the comparison counts pairs.

**Resistance asked of the network, not read from the tissue.** Whether a cell survives depends on where it sits as much as on what it carries. Asking each network how it would respond, at one fixed level, separates the two, and the third route in the table above is exactly the case where they differ.

**All at once.** The Monte Carlo step of one run uses one core (chapter 12); a campaign puts one run on each, which is how forty runs of 3000 steps take about forty minutes on twenty threads.

## Where in the code

`src/campaign.rs`: `campaign`, `parse_seeds`, the outcome of a run and the paired summary. `crates/engine/src/evolution.rs`: `Simulation::relative_responses` and `Simulation::resistant_fraction`, `RESISTANT_BELOW`. `crates/engine/src/network.rs`: `Network::apoptosis_under`, which asks a copy of a network how it would respond. `crates/engine/src/simulation.rs`: `Simulation::pristine`, the unmutated network kept as the yardstick.

`examples/cord-evolving-stress.toml` and `examples/cord-evolving-stress-periodic.toml`, the two arms.

## Still a placeholder

Everything chapters 6 and 8 list, and the measure of resistance: half the source's level, and a half as the line between sensitive and resistant. Twenty seeds give counts, not rates; whether the split is closer to a third or a half needs more. And one schedule of each kind: whether some other rhythm does better is a question the same command can ask, one configuration at a time.
