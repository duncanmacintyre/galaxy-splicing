This directory contains code for computing the protocluster mass budget and preparing Makegalaxy input files.

Key files in this directory:

* *compute_galaxy_mass.ipynb* uses observations to work out the mass budget that we use to guide our SPT2349-56 initial conditions. It outputs galaxy_masses.txt (masses of all galaxies) and galaxy_masses_105.txt (masses of galaxies within 105 kpc of the centre—these are the galaxies used in the simulations).
* *gen_files.py* is a script that generates input files for Makegalaxy from galaxy_masses_105.txt.
* *plot_positions.py* plots a map of the galaxies in the protocluster by reading galaxy_masses.txt and galaxy_masses_105.txt. (This map thus represents observations.)

Note: galaxy_masses_105.txt gives the mass budget according to observations; files like `../splice/chosen-1e6/real010/masses.txt` give the actual mass budget that was used in initial conditions. They may differ slightly. You probably want to use the latter when preparing plots or tables of the simulation's mass budget.