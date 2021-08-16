# plot_snapshot_stellar_mass_histograms.py
# DM: initial version June 2021, plotting part based on Jupyter Notebook by Doug Rennehan
#
#   given GIZMO snapshot files, generate a PDF for each that shows a 2D histogram of the star mass distribution
#
#   for usage run `$ python plot_snapshot_stellar_mass_histograms.py --help`
#   

import multiprocessing
import sys
import os
from argparse import ArgumentParser

import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import h5py

from common import locate_peak_density_3D_and_plot, code_time_to_Myr, grab_property

parser = ArgumentParser(
    description='Given GIZMO snapshot files, generates a PDF for each that shows a 2D histogram of the stellar mass distribution.',
    epilog='The plot for /path/to/some_snap.hdf5 will be saved at out_dir/some_snap.pdf (or with appropriate extension for the format). Required modules: hdf5, scipy-stack.'
)
parser.add_argument('out_dir', help='path to directory in which to save the plots; will be created if it does not already exist')
parser.add_argument('snaps', nargs='+', metavar='snap.hdf5', help='paths to GIZMO snapshot files to plot')
parser.add_argument('--squish-along', nargs='+', dest='squish_along', choices=['x', 'y', 'z'], default=('z',), help='which axis to flatten along when making 2D histograms – use more than one to get multiple subplots; default z')
parser.add_argument('--title', help='title for plots (same for all plots); if not specified, uses snapshot times for titles')
parser.add_argument('--format', default='pdf', help='file save format, e.g. pdf, png, etc.; default pdf')
parser.add_argument('--dpi', default=None, type=int, help='figure dpi (pixels per inch)')
parser.add_argument('--radius', default=75, type=float, help='radius of cube region to plot in kpc; default 75')
parser.add_argument('--nticks', default=7, type=int, help='how many tick marks on plots; default 7')
parser.add_argument('--nbins', default=512, type=int, help='how many histogram bins in each direction; default 512')
args = vars(parser.parse_args())

# directory in which to store PDFs
out_dir = args['out_dir']
# paths to snapshot files
snaps = args['snaps']
# directions to squish along, length gives number of subplots
squish_along = args['squish_along']
# title to use in plots, or None to generate titles from times
title = args['title']
# file save format, e.g. "pdf"
save_format = args['format']
# dpi: integer or None
dpi = args['dpi']
# radius of plot region in kpc
radius = args['radius']
# how many ticks to place on the axes
nticks = args['nticks']
# how many histogram bins in each direction
nbins = args['nbins']

# maps 'x', 'y', 'z' to 0, 1, 2
map_xyz_012 = {'x':0, 'y':1, 'z':2}
# maps name of squished direction to names of directions not squished 
map_not_squished = {'x':('y','z'), 'y':('x','z'), 'z':('x','y')}

# ------ functions

# given a file name for a snapshot, save plot to out_dir
# most of this function's code is by Doug Rennehan
def plot_and_save_from_snapshot_file(fname):

    # load the data from the HDF5 snapshot file
    with h5py.File(fname) as f:
        coords = np.concatenate(
                (grab_property(f, 2, 'Coordinates'),
                 grab_property(f, 3, 'Coordinates'),
                 grab_property(f, 4, 'Coordinates')),
                axis=0)
        masses = np.concatenate(
                (grab_property(f, 2, 'Masses'),
                 grab_property(f, 3, 'Masses'),
                 grab_property(f, 4, 'Masses')),
                axis=0)
        time = code_time_to_Myr(f['/Header'].attrs['Time']) # time in Myr
    
    # set up a figure and axes
    if dpi is not None:
        fig, ax = plt.subplots(1, len(squish_along),
                               figsize=(4*len(squish_along), 3.5),
                               constrained_layout=True,
                               dpi=dpi)
    else:
        fig, ax = plt.subplots(1, len(squish_along),
                               figsize=(4*len(squish_along), 3.5),
                               constrained_layout=True)

    # convert ax to an iterable in the case of one subplot
    if not hasattr(ax, '__iter__'):
        ax = (ax,)

    # plot histogram onto the axes
    locate_peak_density_3D_and_plot(coords, axes=ax, cube_radius=radius, nbins=nbins, 
                                    weights=masses, nticks=nticks,
                                    squish_along=tuple(map_xyz_012[sa] for sa in squish_along))

    # set axis labels
    for a, sa in zip(ax, squish_along):
        axes_names = map_not_squished[sa]
        a.set_xlabel('{} [kpc]'.format(axes_names[0]))
        a.set_ylabel('{} [kpc]'.format(axes_names[1]))

    # set figure title
    fig.suptitle('t = {:.1f} Myr'.format(time) if title is None else title)

    # save figure to disk
    fig.savefig(os.path.join(out_dir, os.path.splitext(os.path.basename(fname))[0] + '.' + save_format))

    # delete figure data from memory to free up space
    fig.clear()
    plt.close(fig)

    print('Finished {}'.format(fname))

# ----- beginning of main script

# create the output directory if it does not already exist
try:
    os.mkdir(out_dir)
except FileExistsError:
    pass

# iterate over all of the snapshots, plot and save each
ncpus = int(os.environ.get('SLURM_CPUS_PER_TASK', default=1))
with multiprocessing.Pool(processes=ncpus) as pool:
    pool.map(plot_and_save_from_snapshot_file, snaps, chunksize=1)
