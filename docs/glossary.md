# Glossary

Every word the course uses, in one place. Chapter numbers point to where each is explained.

**Adhesion** (1). Cells sticking to each other through surface proteins. Modelled as a lower contact energy between cells of the same kind.

**Anoxia** (2). Oxygen so low the cell cannot make energy. Below 15% of the source's level in cellweave, damage to survival begins.

**Apoptosis** (2, 4). Programmed, orderly cell death that a cell can trigger itself. Not modelled.

**Basal layer** (7). The bottom layer of an epithelium, the only one that divides. Cells leave it upwards when they differentiate.

**Bath** (3). A large stirred volume of liquid around a tissue whose oxygen never drops. Modelled by holding every medium pixel the bath can reach at 1.0.

**BCL2** (2). B-cell lymphoma 2, a survival protein. When it falls below 0.25 the cell is necrotic.

**Capillary** (3). The finest blood vessel, where oxygen leaves the blood. No cell in the body is further than 100 to 200 micrometres from one.

**Cell competition** (7). Cells eliminating their neighbours by outdoing them. Here cells with higher ERK press their neighbours small, and the squeezed ones differentiate and leave.

**Cellular Potts model, CPM** (1). Cells as sets of pixels on a lattice, moving pixel by pixel to lower an energy.

**Clonal sweep** (6). One clone outcompeting the rest and taking over the population.

**Clone** (6). Every cell descended from one mutated ancestor with no further mutation in between. The first cell founds clone 0.

**Contact energy, J** (1). What a piece of boundary between two cells of given kinds costs.

**Contact inhibition** (7). A cell with no room around it behaving differently from one with room: here, differentiating more readily.

**Convergent evolution** (6). The same outcome reached through different changes. Seven evolving cords out of eight outgrew their steady state, by at least three routes.

**Cord, tumour cord** (5). Tumour tissue living around a vessel out to the depth oxygen reaches, dead beyond. A steady state with turnover.

**Crowding** (7). How squeezed a cell is: how far below its reference size it sits, as a fraction of it. An input every network can read, like oxygen.

**Determinism** (1, 8). The same seed gives the same simulation, byte for byte. Checked on every push.

**Differentiation** (7). A cell giving up dividing to become a specialised cell. In an epithelium, a basal cell that commits leaves the layer. Here, a chance per step set by the network's differentiation node.

**Diffusion, D** (3). Oxygen wandering from where there is more to where there is less. D is how fast.

**Division** (4). A cell past 1.4 times its lineage's reference volume is cut along a line through its centre of mass; the daughter inherits a copy of the network, or with the mutation rate a copy with one rate changed.

**Drift** (6). A variant spreading or vanishing by chance rather than because of what it does.

**EGF** (2). Epidermal growth factor, the growth signal. A hardcoded input for every cell today.

**Emergence** (6, 7, 8). Behaviour that follows from the model without any rule naming it. The necrotic core, the compact cell, the cord's thickness.

**Energy, Hamiltonian, H** (1). The number every configuration of the lattice has, which the Monte Carlo lowers.

**Epithelium** (7). Tissue made of sheets of cells, lining the skin, the gut, the oesophagus. Where most cancers start.

**Equivalent radius** (5). The radius of a disc with the same area as a shape. How every radius in a tissue's measurements is taken.

**ERK** (2). Extracellular signal-regulated kinase, the last relay of the growth pathway. Sets the target volume.

**Field** (3). One value per pixel, here the oxygen concentration.

**Gone** (4, 7). A cell that has left the tissue for good: resorbed to nothing after dying, or shed after differentiating. Its slot is empty; the lattice keeps a tombstone.

**Growth factor** (2). A molecule other cells release meaning "grow". EGF here.

**HIF-1α** (2). Hypoxia-inducible factor 1 alpha, the protein that accumulates as oxygen falls. Above 0.45 the cell is quiescent.

**Hitchhiking** (6). A neutral mutation spreading because it shares a lineage with one that pays off.

**Hypoxia** (2, 3). Too little oxygen, but enough to live on. Distinct from anoxia.

**Jacobi** (3). The plain way to relax a field, one full update from the previous state. Replaced by SOR for steady states.

**Kinase** (2). An enzyme that activates another protein by attaching a phosphate group.

**Lattice** (1). A regular grid, the word from garden trellis. Here 200 by 200 pixels.

**Lineage reference volume** (4). The size a lineage divides at, inherited unchanged. How a cell type keeps its size across generations.

**Lineage tracing** (7). Marking cells with a heritable label and following what their descendants become. A sheet labels every cell as its own clone.

**Medium** (1, 3). Empty pixels, identity 0. Liquid in a bath; stroma around a vessel.

**Monte Carlo step, MCS** (1). One attempt per pixel to change hands, 40,000 on the default grid. The lattice's unit of time.

**Muller plot** (6). Every clone's share of the population over time, stacked, so a sweep is a band that widens.

**Mutation** (6). A heritable change. Here one rate of a daughter's network scaled by a factor between one half and two.

**Mutation rate** (6). The probability that a daughter is born with a mutation. One in twenty in the evolving cord, a placeholder.

**Necrosis, necrotic** (2, 4). Death by energy failure. A kind on the lattice: stops consuming, network frozen, resorbed.

**Necrotic core** (5). The dead centre of a spheroid too large for oxygen to reach.

**Neutral drift** (7). Clones winning and losing by chance alone, when no cell is better than another. The number of clones falls and the survivors grow.

**NOTCH, NOTCH1** (7). A receptor that pushes a basal cell to differentiate. Its loss lets a clone keep more cells in the basal layer; mutated in much of the normal ageing oesophagus.

**Obstacle** (1, 3). A pixel no cell can take. A vessel is an obstacle held at full oxygen.

**ODE** (2). Ordinary differential equation: rates of change in time only. The network is five of them.

**PDE** (3). Partial differential equation: rates of change in time and space. Oxygen is one.

**Pixel** (1). One square of the lattice, about a micrometre, carrying one identity and one oxygen value.

**Proliferating** (4). Alive, cycling, dividing. The default fate under enough oxygen.

**Quiescent** (4). Alive, not dividing. Hypoxia response above 0.45.

**RAS** (2). The switch at the start of the growth pathway. Stuck on in one cancer in three.

**Red-black** (3). An order for an SOR sweep: the grid coloured like a chessboard, every red pixel updated from its black neighbours, then every black one from the new reds. Each colour can be updated in parallel.

**Relaxation** (3). Iterating a field until it stops moving: its steady state.

**Resorption** (4). Dead tissue being cleared. Two pixels per step here.

**RK4** (2). Fourth order Runge-Kutta, the scheme that advances the network's equations.

**Seed** (1). The number that starts the random draws. Same seed, same run.

**Selection** (6). Variants leaving more descendants because of what they do where they are. Nobody chooses the winner.

**Snapshot** (4, 8). A file listing every living cell's state and the oxygen it reads, one per interval.

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

**Viable depth** (5). How far living tissue reaches from where its oxygen comes in: from the edge inward in a spheroid, from the wall outward around a vessel.

**Volume constraint, λ** (1). The cost of being a different size than the target, times a stiffness lambda.
