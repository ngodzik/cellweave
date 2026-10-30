# Changelog

Versions follow [semantic versioning](https://semver.org). While the version is 0.x the project is in early development: anything may change from one version to the next, configurations and results included. Each version says what changed in what a run does, so that a result can be tied to the version that produced it.

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
