# Changelog

Versions follow [semantic versioning](https://semver.org). While the version is 0.x the project is in early development: anything may change from one version to the next, configurations and results included. Each version says what changed in what a run does, so that a result can be tied to the version that produced it.

## Unreleased

- **A fuller cell.** `examples/networks/fuller-cell.toml` fills in what the RAS/ERK network and the stress response leave out: EGFR, RAF and MEK before ERK; cyclin D, RB, E2F, cyclin E and the brake p27 between hypoxia and the arrest; MDM2, p21 and BAX around p53. Its rates are set so that a cell asks for the same size, stops dividing at the same oxygen, dies of anoxia at the same oxygen and of a stress at the same rate as in the simpler networks, which tests hold (chapter 28). `spheroid-fuller.toml` and `cord-stress-fuller.toml` run it. Under a stress the cord dies out where the simpler one holds, through p21. No other run changes.
- **A term that switches off.** `inhibited_by` joins the catalogue of terms, the mirror of `activated_by`: `- rate * source * self`. No run changes.
- **One cell opened.** In the viewer, a click on a cell opens a window that draws its network as it stands, what it reads, and what it decides; `--open`, `--open-at`, `--follow`, `--protein` and `--window` do the same from the command line, for recordings. The camera turns only while the window has the pointer and the focus.
- **The course.** Chapter 28 is the fuller cell and the window; what the model does not say is now chapter 29.

## 0.2.0 (2026-11-01)

A tissue that grows in a volume, and a spheroid in it that stops near the size measured spheroids reach, once its dead core breaks down and the medium of its well is still. Experimental.

One change is incompatible: a network file that still sets `j_medium`, `j_self_base` or `j_self_from_growth` is refused, since nothing ever read them (first item below).

- **Contact energies are not a network's.** A network file's `[mechanics]` carried `j_medium`, `j_self_base` and `j_self_from_growth`, read into every cell and handed to nothing: the lattice has always priced a boundary by the kinds of the cells across it (chapter 1). They are removed from the example networks, and a network that sets them, or any other unknown key under `[mechanics]`, is now refused with the key named. No run changes.
- **A tissue in a volume.** A configuration's `depth` stacks slices: a tumour in a bath or around a vessel grows in three dimensions, cells dividing across planes, a vessel running through the box as a cylinder, the bath coming in from six faces, radii those of balls (chapter 26). A plane, one slice deep, runs exactly as it did: every example gives the same output, byte for byte. `max_cells` replaces the fixed ceiling of 2000 living cells.
- **Faster.** The field is relaxed and the bath flooded only around the tissue, and in a volume the Monte Carlo step draws only there and the field is over-relaxed at its optimum: the spheroid's first 120 steps on one core take 11.5 seconds instead of 65, with the same result.
- **A dead core that breaks down.** A configuration's `debris_loss` is the share of its volume a dead cell's debris loses each step, as in Greenspan's model of a spheroid: the living tissue closes in over it at the edge of the dead tissue, and gives as much up at its surface. `sweeps_per_step` sweeps the lattice several times a step, so that the tissue carries what the core loses out to its surface. `examples/spheroid-3d.toml` grows the spheroid in a volume with both, and `scripts/compare-spheroids.py` sets it against the measured spheroids as balls, in micrometres, at their end and at their pace (chapter 27). Without `debris_loss` and with one sweep a step, every run is the run it was.
- **A cycle that takes time.** A configuration's `[cycle]` has each cell draw its cycle's duration at birth, log-normal around a median, grow through it from half of what its network asks for to all of it, and divide only once it has passed: a tissue grown in a volume no longer divides all at once every 20 steps, and a daughter no longer asks at birth for twice her volume. Left out, every run is the run it was.
- **Contact to the corners.** A configuration's `contact_order` counts two pixels as touching across edges and corners as well as faces, for the contact energy, each contact weighted so that a boundary along an axis costs what it cost: what a boundary costs then depends far less on its direction, and a tissue in a volume rounds into a ball instead of squaring into a cube (chapter 27). Copies are still proposed across faces alone. Left at 1, every run is the run it was.
- **Fluid that drains.** With a debris loss, medium the tissue closes in on every side, out of the bath's reach, is closed in by the living tissue around it, one layer a step, as the debris is: a core no longer hollows into a cavity that takes what it loses.
- **A waste that arrests.** A configuration's `[waste]` has every living cell release a waste, Greenspan's inhibitor, that diffuses as oxygen does and that the bath washes away; each cell's network reads it as the input `waste`. `examples/networks/ras-erk-waste.toml` gives the arrest to a node p27 raised by HIF-1α and by the waste (chapter 27). Left out, every run is the run it was; a gained link still reads only the four earlier inputs.
- **A factor the dead release.** A configuration's `[necrotic_factor]` has every dead pixel release a factor that diffuses and that the bath washes away, as Freyer found the necrosis of spheroids does, and each cell's network reads it as `necrotic_factor`. A new term, `damaged_above`, damages a node past a level of its source; `examples/networks/ras-erk-necrotic.toml` has the factor kill through BCL2. Left out, every run is the run it was.
- **A still medium.** A configuration's `[still_medium]` has the bath's medium still, oxygen diffusing through it at its `diffusion`: the level next to a spheroid in a volume is then that of a still medium around a ball, lower the more the spheroid takes up, as measured at the surface of spheroids in culture (chapter 27). The spheroid in a volume runs in one, its diffusion set against the measured spheroids' end, and stops at 368 µm against their 348, its cores growing with size as theirs do, by half as much. Left out, the medium is held at the source's level up to the tissue, and every run is the run it was.
- **A volume in the viewer.** A tissue grown in a volume is drawn as its voxels, cut open through the middle so that its core shows; C opens and closes the cut (chapter 26).
- **Legends on the animations.** Every animation in the README says in its corner what its colours mean, and the README says they are sped up.
- **The spheroid in a volume, in the README.** `examples/spheroid-3d.toml`, recorded in the viewer and cut open, from one ball of cells to its end.

## 0.1.0 (2026-10-30)

The first release: everything built so far, and the first calibration on public data. Experimental.

**The model.** Cells as sets of pixels on a lattice (a Cellular Potts model); inside each, a network of proteins read from a TOML file, integrated cell by cell; between them, oxygen, a lethal stress, VEGF and a chemokine diffusing and consumed where the cells are. Cells grow, arrest, divide, die of anoxia, start apoptosis, or leave a renewing layer, each because its network says so.

**Scenarios**, each a file in `examples/`:

- a spheroid in a bath, starving from the inside; a tumour cord around a vessel; a tumour against a capillary lying in the plane;
- mutation and rewiring of the network, and selection; a lethal stress, continuous or periodic, and the resistance it selects;
- a renewing epithelium whose clones drift and, with mutation, sweep; an early tumour placed in it, purged or not by its neighbours;
- vessels that grow towards the cells that call for them, carry blood once joined or as pressure drives it (Poiseuille's law), deliver oxygen in proportion to their flow, regress and are pruned;
- killer immune cells that arrive through the vessels, climb a chemokine and kill, and a shield (PD-L1) tumour cells can mutate.

**Calibrated on public data**, the data fetched by the user and never stored in the repository (`datasets/`, `scripts/fetch-dataset.py`):

- the tumour's oxygen uptake, on the radius at which 994 melanoma spheroids form a core (Browning et al. 2021), checked on tumour cords (Tannock 1968): every tumour configuration and the program's default now use it;
- immune cells' killing, on measured rates (Weigelin et al. 2021, Halle et al. 2016);
- the renewing layer's NOTCH, the setting nearest the clones counted in the mouse oesophagus (Piedrafita et al. 2020) of those tried; still several times too fast after a month;
- a time unit for the sheet (0.56 hours a step) and the tumour (0.3 to 0.9 hours).

**Tools.** Campaigns over many seeds, paired against a second configuration (`--seeds`, `--versus`); rewinding a history by suppressing or undoing a mutation (`--attribute`); lineage tracing (`[tracing]`); scripts that set runs against measured clones and spheroids.

**The viewer**, a Bevy window over the same simulation: the tissue in relief, coloured by fate, protein, oxygen, stress, VEGF, clone or blood flow; two runs side by side; `--record` saves every step as an image, and `--zoom`, `--spin`, `--no-panel` and `--batch` frame it, to make an animation like those in the README.

**The course**, 26 chapters in `docs/` that follow the code from the biology up, with the predictions written before each experiment and the ones that failed.
