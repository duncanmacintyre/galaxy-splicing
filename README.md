This repository contains code for creating [Gizmo](http://www.tapir.caltech.edu/~phopkins/Site/GIZMO_files/gizmo_documentation.html) initial conditions of galaxy protoclusters. These initial conditions can then be used to run simulations similar to those by [Rennehan+2020](https://arxiv.org/abs/1907.00977v2).

## Overview

There are three main steps in making such initial conditions:

1. **Mass budget.** We must decide how many galaxies our protocluster will include and choose a stellar mass, gas mass, and dark matter mass for each. We must also choose our simulation resolution in this step.<sup>[1](#footnote1)</sup>
2. **Create initial conditions for each galaxy.** We use the Makegalaxy software to generate initial conditions for each galaxy.<sup>[2](#footnote2)</sup> Then, we use Gizmo to simulate each galaxy in isolation to ensure that its structure is stable.<sup>[3](#footnote3)</sup>
3. **Create initial conditions for the protocluster.** Once we have ICs for each individual galaxy, we splice them together to form initial conditions for the protocluster.

The code is currently set up to form initial conditions for objects similar to SPT2349-56. We make the mass budget as similar to SPT2349-56 as possible, but randomize the initial position and velocity for each galaxy. Different combinations of these galaxy positions and velocities yield different "realizations".

We could perform steps 1 and 2 once but repeat step 3 multiple times to create several realizations of the protocluster while keeping the same mass budget. Conversely, we might find it useful to start with various mass budgets but always keep the same galaxy positions and velocities to see how a realization's evolution is affected by the mass budget.

## Directories

The *generate* directory contains calculations that work out the mass budget from observations (compute_galaxy_mass.ipynb). It also has code for creating Makegalaxy input files (gen_files.py). We store our Makegalaxy output in the makegalaxy-1e6 and makegalaxy-1e7 subdirectories.

The *splice* directory contains code to easily start Gizmo simulations of galaxies in isolation; we store the results for each galaxy in the chosen-1e6 and chosen-1e7 subdirectories. The generate_protocluster.py script splices these results together to form initial conditions for the protocluster. Results for several realizations of splicing are included, for example, at `splice/chosen-1e7/real010.hdf5`.

The *analysis* directory contains some scripts for analyzing simulation results. They are also useful for inspecting initial conditions files as they are being made. For instance, the plot_snapshot_stellar_mass_histograms.py is useful for visualizing the galaxies as they are simulated in isolation.


---

<a name="footnote1">1</a>: "Resolution" just means the initial mass of each particle. If we set each gas particle to have an initial mass of 1e6 M☉, we will have ten times more particles and a better resolution than if we had set each gas particle to have an initial mass of 1e7 M☉. We make dark matter particles five times more massive than other types because of the large amount of dark matter and because we are more interested in other particle types. Thus, I might say "1e6 resolution" to mean that gas particles have initial masses of 1e6 M☉, start particles have initial masses of 1e6 M☉, and dark matter particles have initial masses of 5e6 M☉.

<a name="footnote2">2</a>: The Makegalaxy software is not included in this repository.

<a name="footnote3">3</a>: We instruct Makegalaxy to form initial conditions with gas fractions *higher* than those we desire. Then, we simulate in isolation until the desired gas fraction is reached (as gas converts to stars due to star formation). You can see how the boost is applied in the generate/generate_makegalaxy_params.py file. After simulating in isolation, we must choose which output snapshot gives us the gas fraction closest to what we desire; this is done by the splice/choose-snaps.py script. The chosen snapshots are in the splice/chosen* directories.
