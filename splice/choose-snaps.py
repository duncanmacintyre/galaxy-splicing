# choose snapshots for splicing
# call this script with arguments of simulation directory names; requres modules: hdf5, scipy-stack
#
# This script iterates over each directory and
#   1. Chooses which snapshot we should use for the initial conditions;*
#   2. Copies the .hdf5 file to an output directory;
#   3. Plots the star formation rate, with an arrow indicating the chosen snapshot, and stats; and
#   4. Saves a PDF of this plot to the output directory.
#
# *strategy: we choose first snapshot satisfying either
#          (1) less than 10% change in SFR over subsequent 10 Myr and less than 3.33% change in
#              SFR over subsequent 5 Myr, or
#          (2) SFR remains below 2% of the total peak value for subsequent 10 Myr
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


# ----- constants -----

code_mass_units_in_Msun = 1e10
code_time_units_in_Myr = 978.028

time_between_snaps = 0.002 # time between snapshots in code units
time_snapshot_zero = 0 # time of snapshot 0 in code units

out_dir = 'chosen' # directory to save chosen .hdf5 files and SFR plots in

# strategy: we choose first snapshot satisfying either
# (1) less than 3.33% and 10% change in SFR over subsequent 5 Myr and 10 Myr respectively, or
# (2) SFR is less than 2% of the peak value for subsequent 10 Myr
# interval = 10 # size of interval in which we consider change, in Myr
cutoff = 0.1 # fractional change must be less than this amount
frac_accepted_below = 0.02
# SFR must change less than cutoff over interval and less than cutoff/3 over interval/2,
# or must be less than frac_accepted_below for at least interval

# ----- functions -----

# given non-negative time in Myr, return the number of first snapshot that would be at or after that time
def time_to_snap(t):
    return math.ceil((t/code_time_units_in_Myr - time_snapshot_zero) / time_between_snaps)
assert time_to_snap(time_snapshot_zero*code_time_units_in_Myr) == 0
assert time_to_snap((time_snapshot_zero + 0.01*time_between_snaps)*code_time_units_in_Myr) == 1
assert time_to_snap((time_snapshot_zero + 6*time_between_snaps)*code_time_units_in_Myr) == 6
assert time_to_snap(20) == math.ceil((20/code_time_units_in_Myr - time_snapshot_zero) / time_between_snaps)

# given a snapshot number, return its time in Myr
def snap_to_time(n):
    return (time_snapshot_zero + n * time_between_snaps) * code_time_units_in_Myr
assert snap_to_time(20) == (time_snapshot_zero + 20 * time_between_snaps) * code_time_units_in_Myr

# get time and SFR from file of given path sfr.txt in working directory and return as tuple of two arrays in Myr, Msun/yr respectively
def load_sfr(sfr_file_path):
    df = pandas.read_csv(sfr_file_path, sep=' ', names=('t', 'expected', 'SFR', 'M* per step', 'M* this step'))
    t = df['t'].to_numpy() * code_time_units_in_Myr
    sfr = df['SFR'].to_numpy() * code_mass_units_in_Msun / (code_time_units_in_Myr * 1e6)
    return (t, sfr)

# given time, sfr, snapshot's path and time in Myr, initial conditions path, and a galaxy name: print stats for the snapshot, plot the SFR with stats, and save figure to given path
def plot_sfr(t, sfr, snapshot_path, snapshot_time, ic_path, galaxy_name, figure_save_path):
    snapshot_name = os.path.basename(snapshot_path)
    desc = make_description_of_gas_lost(snapshot_path, ic_path)
    print('{} {}: {}'.format(galaxy_name, snapshot_name, desc))
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

# given paths to snapshot and initial conditions: return string with how much gas lost to star formation
def make_description_of_gas_lost(snapshot_path, ic_path):
    N_snap, m_snap = get_amount_gas_particles(snapshot_path)
    N_zero, m_zero = get_amount_gas_particles(ic_path)
    return '{:d} of {:d} gas particles lost ({:.1f}%, {:.1f}% by mass)'.format(
                N_zero-N_snap, N_zero, 100*(N_zero-N_snap)/N_zero, 100*(m_zero-m_snap)/m_zero)

# given path to snapshot, return its number of gas particles and total gas mass in code units
def get_amount_gas_particles(snapshot_path):
    with h5py.File(snapshot_path) as f:
        masses = f['/PartType0/Masses']
        return (len(masses), np.sum(masses))

# given non-negative int of snapshot number, return filename of its HDF5 file (with extension, without full path)
def get_snap_filename(n):
    return 'snapshot_{:0>3d}.hdf5'.format(n)
assert get_snap_filename(0) == 'snapshot_000.hdf5'
assert get_snap_filename(25) == 'snapshot_025.hdf5'

# given non-negative int of snapshot number, return full path to its HDF5 file in the snap_dir folder
def get_snap_path(n, snap_dir):
    return os.path.join(snap_dir, get_snap_filename(n))

# given time and SFR arrays in Myr and Msun/yr, return number for the snapshot we should choose for initial conditions
def choose_snap(t, sfr):
    # find time when galaxy is considered "stable" as per our strategy
    # see also https://stackoverflow.com/a/8534381/13326516
    accept_below = frac_accepted_below*np.max(sfr)
    interval = np.max(t)/20
    chosen_time = next(
            (this_t for this_t, this_sfr in zip(t, sfr) if
                within_cutoff(sfr[(t > this_t) & (t < this_t+interval) & (sfr != 0)],
                              this_sfr, cutoff, accept_below)
                and within_cutoff(sfr[(t > this_t) & (t < this_t+interval/2) & (sfr != 0)],
                                  this_sfr, cutoff/3, accept_below)
                and within_cutoff(sfr[(t > this_t) & (t < this_t+interval*3) & (sfr != 0)],
                                  this_sfr, cutoff*2, accept_below)
                and within_cutoff(sfr[(t > this_t) & (t < this_t+interval*4) & (sfr != 0)],
                                  this_sfr, cutoff*3.5, accept_below)),
            None) # default: will be None if condition never reached
    if chosen_time is not None: # if we found a time, make sure it is not later than last snapshot
        chosen_time = min(chosen_time, t[-1] - time_between_snaps*code_time_units_in_Myr)
    else: # if we could not find a time, default to the last snapshot
        chosen_time = t[-1] - time_between_snaps*code_time_units_in_Myr
    return time_to_snap(chosen_time)
# Note: This method of making sure the snapshot number does not exceed the last snapshot is probably
# a bit buggy. It would be better to make a separate function to figure out what the max snapshot 
# number is by looking at what files exist.
# Practically, this doesn't matter much because we expect to choose a time long before the last snapshot.

# return True if either
#  (1) test_val>0, len(x)>0, and all values of x within cutoff fraction of test_val, or
#  (2) test_val>0, len(x)>0, and all values of x have magnitude less than accept_below;
# otherwise return False
def within_cutoff(x, test_val, cutoff, accept_below=0):
    return (test_val > 0 and len(x) > 0
            and (np.all(x/test_val < (1 + cutoff)) and np.all(x/test_val > 1/(1 + cutoff)))
                or np.all(np.abs(x) < accept_below))

# given path to galaxy sim directory, choose snapshot, plot, and save output
def process_one_sim(this_dir):
    t, sfr = load_sfr(os.path.join(this_dir, 'sfr.txt'))
    n = choose_snap(t, sfr)
    plot_sfr(t=t, sfr=sfr, snapshot_path=get_snap_path(n, this_dir), snapshot_time=snap_to_time(n),
             ic_path=get_snap_path(0, this_dir), galaxy_name=this_dir,
             figure_save_path=os.path.join(out_dir, this_dir + '.sfr.pdf'))
    shutil.copyfile(get_snap_path(n, this_dir), os.path.join(out_dir, this_dir + '.hdf5'))



# ----- script part -----

# iterate over each argument given, running with multiple cores if they are given with slurm, one core otherwise
ncpus = int(os.environ.get('SLURM_CPUS_PER_TASK', default=1))
print('Starting to choose snapshots with {} CPUs...'.format(ncpus))
with multiprocessing.Pool(processes=ncpus) as pool:
    pool.map(process_one_sim, sys.argv[1:])
    

