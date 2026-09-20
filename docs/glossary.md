# Glossary

Every word the course uses, in one place. Chapter numbers point to where each is explained.

**Adhesion** (1). Cells sticking to each other through surface proteins. Modelled as a lower contact energy between cells of the same kind.

**Anoxia** (2). Oxygen so low the cell cannot make energy. Below 15% of the source's level in cellweave, damage to survival begins.

**Apoptosis** (2, 4). Programmed, orderly cell death that a cell can trigger itself. Not modelled.

**Bath** (3). A large stirred volume of liquid around a tissue whose oxygen never drops. Modelled by holding every medium pixel the bath can reach at 1.0.

**BCL2** (2). B-cell lymphoma 2, a survival protein. When it falls below 0.25 the cell is necrotic.

**Capillary** (3). The finest blood vessel, where oxygen leaves the blood. No cell in the body is further than 100 to 200 micrometres from one.

**Cellular Potts model, CPM** (1). Cells as sets of pixels on a lattice, moving pixel by pixel to lower an energy.

**Contact energy, J** (1). What a piece of boundary between two cells of given kinds costs.

**Cord, tumour cord** (5). Tumour tissue living around a vessel out to the depth oxygen reaches, dead beyond. A steady state with turnover.

**Determinism** (1, 6). The same seed gives the same simulation, byte for byte. Checked on every push.

**Diffusion, D** (3). Oxygen wandering from where there is more to where there is less. D is how fast.

**Division** (4). A cell past 1.4 times its lineage's reference volume is cut along a line through its centre of mass; the daughter inherits a copy of the network.

**EGF** (2). Epidermal growth factor, the growth signal. A hardcoded input for every cell today.

**Emergence** (6). Behaviour that follows from the model without any rule naming it. The necrotic core, the compact cell, the cord's thickness.

**Energy, Hamiltonian, H** (1). The number every configuration of the lattice has, which the Monte Carlo lowers.

**ERK** (2). Extracellular signal-regulated kinase, the last relay of the growth pathway. Sets the target volume.

**Field** (3). One value per pixel, here the oxygen concentration.

**Gone** (4). A necrotic cell that has been resorbed to nothing. Its slot is empty; the lattice keeps a tombstone.

**Growth factor** (2). A molecule other cells release meaning "grow". EGF here.

**HIF-1α** (2). Hypoxia-inducible factor 1 alpha, the protein that accumulates as oxygen falls. Above 0.45 the cell is quiescent.

**Hypoxia** (2, 3). Too little oxygen, but enough to live on. Distinct from anoxia.

**Jacobi** (3). The plain way to relax a field, one full update from the previous state. Replaced by SOR for steady states.

**Kinase** (2). An enzyme that activates another protein by attaching a phosphate group.

**Lattice** (1). A regular grid, the word from garden trellis. Here 200 by 200 pixels.

**Lineage reference volume** (4). The size a lineage divides at, inherited unchanged. How a cell type keeps its size across generations.

**Medium** (1, 3). Empty pixels, identity 0. Liquid in a bath; stroma around a vessel.

**Monte Carlo step, MCS** (1). One attempt per pixel to change hands, 40,000 on the default grid. The lattice's unit of time.

**Necrosis, necrotic** (2, 4). Death by energy failure. A kind on the lattice: stops consuming, network frozen, resorbed.

**Necrotic core** (5). The dead centre of a spheroid too large for oxygen to reach.

**Obstacle** (1, 3). A pixel no cell can take. A vessel is an obstacle held at full oxygen.

**ODE** (2). Ordinary differential equation: rates of change in time only. The network is five of them.

**PDE** (3). Partial differential equation: rates of change in time and space. Oxygen is one.

**Pixel** (1). One square of the lattice, about a micrometre, carrying one identity and one oxygen value.

**Proliferating** (4). Alive, cycling, dividing. The default fate under enough oxygen.

**Quiescent** (4). Alive, not dividing. Hypoxia response above 0.45.

**RAS** (2). The switch at the start of the growth pathway. Stuck on in one cancer in three.

**Relaxation** (3). Iterating a field until it stops moving: its steady state.

**Resorption** (4). Dead tissue being cleared. Two pixels per step here.

**RK4** (2). Fourth order Runge-Kutta, the scheme that advances the network's equations.

**Seed** (1). The number that starts the random draws. Same seed, same run.

**Snapshot** (4, 6). A file listing every living cell's state and the oxygen it reads, one per interval.

**SOR** (3). Successive over-relaxation, the sweep that solves for a steady state in about L passes rather than L².

**Spheroid** (5). A ball of tumour cells grown in liquid. Starves from the inside past a few hundred micrometres.

**Steady state, equilibrium** (3). A field that has stopped changing. Where oxygen always is, at the scale of a cell's decisions.

**Stroma** (3). Ordinary tissue between and around tumour cells. Consumes too.

**Target volume** (1, 4). The size a cell tries to be, set by its ERK.

**Temperature, T** (1). How often a costly move is accepted anyway. The noise knob.

**Turnover** (5). Births balancing deaths in a population of steady size.

**Uptake, u** (3). What a pixel of tissue takes up per unit time. A sink in the field's equation.

**VEGF** (2). Vascular endothelial growth factor, the call for new vessels. Computed, unused.

**Vessel** (3, 5). An obstacle disc held at full oxygen in the middle of the box.

**Volume constraint, λ** (1). The cost of being a different size than the target, times a stiffness lambda.
