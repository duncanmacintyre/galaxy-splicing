This directory contains scripts for data analysis.

Key scripts:

### common.py

This file defines some functions used in other files. Of particular interest is the common.grab_property() function. It is the best way to extract data from HDF5 files (e.g. snapshots, initial conditions) while properly handling missing data. For example, the following loads the coordinates for PartType4 (stars):

```
>>> import h5py
>>> import common
>>> f = h5py.File('some_snapshot_or_IC.hdf5', 'r')
>>> gas_mass = common.grab_property(f, 4, 'Coordinates')
```

If there are no particles of PartType4, this will return an empty array of the correct shape rather than giving an error.

### plot_snapshot_stellar_mass_histograms.py

Run `$ python script_name.py --help` to see usage.

This script plots 2D histograms of the stellar mass from Gizmo snapshot files (or other HDF5 files, e.g. initial conditions files). These plots are like maps. One useful technique is to run this script on all snapshots for a certain simulation then combine the resulting plots together to form a rudimentary movie, with each snapshot giving one frame. (One way to combine the plots on macOS is with QuickTime's [image sequence feature](https://support.apple.com/en-ca/guide/quicktime-player/qtp315cce984/mac).) Example usage:

```
$ python plot_snapshot_stellar_mass_histograms.py ./ snap* --format png --dpi 200 --radius 120 --squish-along y z x
```

![example output of plot_snapshot_stellar_mass_histograms.py](../splice/chosen-1e6/real010-log/real010_no_arrows.pdf)

This script is multi-threaded; you can use Slurm scripts like:

```
#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --mem=32000
#SBATCH --cpus-per-task=8
module load hdf5 scipy-stack
python /path/to/plot_snapshot_stellar_mass_histograms.py ./ snap* --format pdf --dpi 200 --radius 120 --squish-along y z x
```

### get_masses_from_snapshots.py

Run `$ python script_name.py --help` to see usage.

This script goes through snapshots and tabulates the mass contained within radial bins. This results in several mass table .txt files (unless you use the `--print` option). You can see example output at `analysis/stellar_mass_radial_envelope/SPT2349_1e6_gf0.9/` (the mass_table* files).

This script is multithreaded; I like to use Slurm scripts like:


```
#!/bin/bash
#SBATCH --time=00:30:00
#SBATCH --mem=32000
#SBATCH --cpus-per-task=8
module load hdf5 scipy-stack
python /path/to/get_masses_from_snapshots.py -f snapshot_{000...400} -R $(cat /path/to/radius_bins.txt) 
```

This will operate on files snapshot_000.hdf5, snapshot_001.hdf5, ..., snapshot_400.hdf5 in the directory in which it is submitted. Observe the radius_bins.txt file; I often `cat` it to get consistent radii to use with the `-R` argument. (The Gizmo parameter file defines the time between snapshots in code units.)

### stellar_mass_radial_envelope/plot_mass_envelopes.py

This file defines a function which makes nice plots of how the radial distribution of mass changes over time. Example output can be seen at `analysis/stellar_mass_radial_envelope/SPT2349_1e6_gf0.9/` (the PDF files).

Note that this is not a script; the function must be imported and run elsewhere. An example of how to do this is at `analysis/stellar_mass_radial_envelope/SPT2349_1e6_gf0.9/plot_Mstar_envelope.py`. This script takes the output tables from get_masses_from_snapshots.py and passes them to plot_mass_envelopes.py.

![example output of plot_mass_envelopes.py](stellar_mass_radial_envelope/SPT2349_1e6_gf0.9/org_fgas=0.9_mean_absolute_Mstar_envelope.pdf)