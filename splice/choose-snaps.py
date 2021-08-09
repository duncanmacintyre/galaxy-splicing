# choose snapshots for splicing
# call this script with arguments of simulation directory names; requres modules: hdf5, scipy-stack
#
# This script iterates over each directory and
#   1. Chooses which snapshot we should use for the initial conditions;*
#   2. Copies the .hdf5 file to an output directory;
#   3. Plots the star formation rate, with an arrow indicating the chosen snapshot, and stats; and
#   4. Saves a PDF of this plot to the output directory.
#
# *strategy: We choose the first snapshot satisfying either
#          (1) SFR remains below 1% of the global maximum for the 30 Myr after the snapshot.
#    or    (2) Over subsequent 10 Myr, SFR max/min < 10%; over subsequent 15 Myr, SFR max/min < 20%;
#              and over subsequent 45 Myr, SFR max/min < 55%.
#

import sys
import os
import shutil
import math
import multiprocessing

import h5py

import numpy as np
import pandas
import matplotlib.pyplot as plt
import matplotlib.patches

import common


# ----- constants -----

code_mass_units_in_Msun = 1e10
code_time_units_in_Myr = 978.028

time_between_snaps = 0.002 # time between snapshots in code units
time_snapshot_zero = 0 # time of snapshot 0 in code units

out_dir = 'chosen' # directory to save chosen .hdf5 files and SFR plots in

# path to fixed-with file of with desired gas fractions
# must have a column 'label' with galaxy names and a column 'fgas' with gas fractions
path_to_gas_fraction_table = os.path.join(
    os.path.dirname(os.path.abspath(__file__)), '..', 'generate', 'galaxy_masses_105.txt')


# ----- setup -----

# the dictionary fgas_target_dict has keys of galaxy names and values of desired gas fractions
df_fgas = pandas.read_fwf(path_to_gas_fraction_table)
fgas_target_dict = {k:v for k,v in zip(df_fgas['label'], df_fgas['fgas'])}
del df_fgas


# ----- functions -----

# given path to galaxy's Gizmo sim directory, choose snapshot, plot, and save output
def process_one_sim(this_dir):
    try:
        # if this_dir = '/path/to/some_direc', galaxy_name = 'some_direc'
        galaxy_name = os.path.basename(this_dir)
        # gas fraction that we want this galaxy to have
        fgas_target = fgas_target_dict[galaxy_name]
        # choose snapshot that gives desired gas fraction
        #       chosen_snap     its path
        #       fgas            its gas fraction
        #       chosen_time     its time
        chosen_snap, fgas, chosen_time = choose_snap(this_dir, fgas_target)

        # load values for star formation rate and corresponding times
        t, sfr = load_sfr(os.path.join(this_dir, 'sfr.txt'))
        # plot the SFR with an arrow pointing to the chosen snapshot
        plot_sfr(t=t, sfr=sfr,
                 snapshot_path=chosen_snap, snapshot_time=common.code_time_to_Myr(chosen_time),
                 galaxy_name=galaxy_name, fgas=fgas, fgas_target=fgas_target,
                 figure_save_path=os.path.join(out_dir, this_dir + '.sfr.pdf'))

        # copy the snapshot to the output directory
        shutil.copyfile(chosen_snap, os.path.join(out_dir, galaxy_name + '.hdf5'))
    except Exception as e:
        print('An error was encountered while processing {}. Skipping.'.format(this_dir))
        print(repr(e))

# given path to directory with Gizmo snapshots, choose gas fraction closest to fgas_target
#   this_dir is the path to a directory containing Gizmo snapshot files
#   returns tuple with three items:
#       * path to chosen snapshot file inside this_dir, e.g. 'snapshot_023.hdf5'
#       * actual gas fraction of the chosen snapshot file, e.g. 0.72
#       * time of chosen snapshot in code units
def choose_snap(this_dir, fgas_target):
    # get snapshot file names and gas fractions, and throw an error if no snapshot files are found
    snap_paths = find_snapshot_files(this_dir)
    if len(snap_paths)==0:
        print('No snapshot files were found!')
        assert False
    snap_fgas = np.fromiter((get_snapshot_fgas(s) for s in snap_paths),
                            count=len(snap_paths), dtype='float64')
    # get index of the snapshot in snap_paths that has a gas fraction closest to the target
    n = np.argmin(np.abs(snap_fgas - fgas_target))
    # get time of chosen snapshot
    with h5py.File(snap_paths[n], 'r') as f:
        chosen_time = f['Header'].attrs['Time']
    return (snap_paths[n], snap_fgas[n], chosen_time)

# return tuple with paths to snapshot files inside this_dir
def find_snapshot_files(this_dir):
    return tuple((os.path.join(this_dir, s) for s in os.listdir(this_dir) if is_snapshot(s)))

# given a file name, return True if it is a snapshot, False otherwise
# e.g. is_snapshot('some_dir/snapshot_023.hdf5') gives True, is_snapshot('foo') gives False
def is_snapshot(s):
    b = os.path.basename(s)
    return b.startswith('snapshot') and b.endswith('.hdf5')

# return the gas fraction of a snapshot given its path
def get_snapshot_fgas(snapshot_path):
    with h5py.File(snapshot_path, 'r') as f:
        Mgas = np.sum(common.grab_property(f, 0, 'Masses'))
        Mstars = (np.sum(common.grab_property(f, 2, 'Masses'))
                  + np.sum(common.grab_property(f, 3, 'Masses'))
                  + np.sum(common.grab_property(f, 4, 'Masses')))
    return Mgas/(Mgas+Mstars)

# get time and SFR from file of given path sfr.txt in working directory and return as tuple of two arrays in Myr, Msun/yr respectively
def load_sfr(sfr_file_path):
    df = pandas.read_csv(sfr_file_path, sep=' ', names=('t', 'expected', 'SFR', 'M* per step', 'M* this step'))
    t = df['t'].to_numpy() * code_time_units_in_Myr
    sfr = df['SFR'].to_numpy() * code_mass_units_in_Msun / (code_time_units_in_Myr * 1e6)
    return (t, sfr)

# given time, sfr, snapshot's path and time in Myr, initial conditions path, and a galaxy name: print stats for the snapshot, plot the SFR with stats, and save figure to given path
def plot_sfr(t, sfr, snapshot_path, snapshot_time, galaxy_name,
             fgas, fgas_target, figure_save_path):
    snapshot_name = os.path.basename(snapshot_path)
    desc = 'fgas: {:.3f} target: {:.3f}'.format(fgas, fgas_target)
    print('{} {} ({})'.format(galaxy_name, snapshot_name, desc))
    sfr_for_snapshot_time = np.interp(snapshot_time, t, sfr)

    fig, ax = plt.subplots(constrained_layout=True)
    ax.plot(t, sfr)
    ax.annotate('',
                xy=(snapshot_time, sfr_for_snapshot_time),
                xytext=(snapshot_time, sfr_for_snapshot_time+0.1*np.max(sfr)),
                arrowprops={'arrowstyle':'->', 'color':'red'})
    ax.set_xlabel('Time [Myr]')
    ax.set_ylabel('SFR [M☉/yr]')
    ax.set_title('Star formation rate for {}'.format(galaxy_name))

    # display stats using legend
    # code borrowed from https://stackoverflow.com/a/59109053/13326516
    handle = matplotlib.patches.Rectangle((0, 0), 1, 1, fc="white", ec="white", lw=0, alpha=0)
    label = '{}: {}'.format(snapshot_name, desc)
    ax.legend((handle,), (label,), loc='best', fontsize='small', 
              fancybox=True, framealpha=0.7, 
              handlelength=0, handletextpad=0)

    fig.savefig(figure_save_path)
    plt.close(fig)


# ----- script part -----

# iterate over each argument given, running with multiple cores if they are given with slurm, one core otherwise
ncpus = int(os.environ.get('SLURM_CPUS_PER_TASK', default=1))
print('Starting to choose snapshots with {} CPUs...'.format(ncpus))
with multiprocessing.Pool(processes=ncpus) as pool:
    pool.map(process_one_sim, sys.argv[1:])
    


