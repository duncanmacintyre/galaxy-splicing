This directory contains code to (1) help simulate galaxies in isolation and (2) splice the results together to form protocluster initial conditions.

Key scripts:
* *start_gizmo_isolated_galaxy.py* operates on a Makegalaxy HDF5 output file to start simulating the galaxy in isolation. It estimates the resources required, creates a new simulation directory, and starts Gizmo with the Slurm scheduler.
* *choose-snaps.py* chooses HDF5 snapshot files from isolated simulation output that we will use for the final initial conditions for each galaxy.
* *generate_protocluster.py* splices together HDF5 files for different galaxies to form initial conditions for the protocluster. It can assign random positions and velocities for the different galaxies or load already-chosen positions and velocities with the `--load` argument.

Run `$ python script_name.py --help` to see usage.

Other interesting files:

**chosen-1e6/\*.hdf5, chosen-1e7/\*.hdf5** (except for realizations, below)

These are the final initial conditions for each galaxy in 1e6 and 1e7 resolution respectively. They were chosen from the results of simulation in isolation.

**chosen-1e6/real\*.hdf5, chosen-1e7/real\*.hdf5**

These HDF5 files give protocluster initial conditions for certain sets of initial positions and velocities for the galaxies. `chosen-1e6/realXXX.hdf5` and `chosen-1e7/realXXX.hdf5` are the same realization (set of positions and velocities) but with different resolutions. 

**chosen-1e6/real\*-log/, chosen-1e6/real\*-log/**

These directories contain supporting "log" material for the protocluster initial conditions. Each has a table of galaxy masses, PDFs with plots, as well as Gizmo parameter files (.param) and Slurm batch scripts (.slurm_job.sh) with some fields filled in (e.g. masses). These folders also contain pickle files with the actual positions and velocities assigned to the galaxies, as well as .hdf5 symlinks—together they can be used by generate_protocluster.py to reconstruct the initial conditions with the `--load` option.

