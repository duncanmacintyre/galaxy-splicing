# generate_protocluster.py
#
# DR: original version for 2020 paper
# DM: rewrote in June 2021
#
# Splice together individual galaxies to form a cluster. Combines .hdf5 GIZMO
# snapshot files to generate one .hdf5 initial conditions file. For usage, do
#
#       $ python generate_protocluster.py --help
#
# required modules: hdf5, scipy-stack
#
# If you like, you can wrap this inside a shell script like:
#
#       #!/bin/sh
#       module load hdf5 scipy-stack
#       python ~/path/to/generate_protocluster.py "$@"
#

import os
import shutil
import pickle
import argparse
import h5py
import numpy as np
from pandas import DataFrame
import matplotlib.pyplot as plt

import common
from generate_protocluster_functions import *

# ----- constants -----

double_precision = True  # use True if GIZMO expects double precision ICs, False if single precision
do_bulge = False         # whether to include bulge particles in output (we always include disk particles)
do_star = True           # whether to include formed star particles in output (we always include disk particles)
density_present = False  # whether to include the gas field Density in output

# these values for r_maximum, v_std, and v_maximum are
# derived in the file ../generate/position_velocity_distribution.ipynb

# paths to Slurm script template and GIZMO parameter file template
# os.path.dirname(__file__) is the folder containing this script
path_to_slurm_script_template = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                             'protocluster.slurm_job.sh')
path_to_gizmo_param_template = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                            'protocluster-1e7.param')

# ----- setup -----

parser = argparse.ArgumentParser(description='Splice together individual galaxies to form a cluster.\n Combines .hdf5 GIZMO snapshot files to generate one .hdf5 initial conditions file.',)
parser.add_argument('output_file', type=str, help='path to output HDF5 file to generate')
parser.add_argument('input_file', nargs='*', type=str, help='paths to HDF5 files to splice together')
parser.add_argument('--log', metavar='LOG_DIR', help='directory in which to save log files; if not specified we use output_file with .log added to the end')
parser.add_argument('--seed', type=int, default=None, help='integer seed for generating random variables; omit for random seed')
parser.add_argument('--no-copy', action='store_true', dest='no_copy', help="don't copy the input files into the log folder, use symbolic links instead")
parser.add_argument('--load', metavar='LOAD_DIR', help='load a previously stored configuration; LOAD_DIR is the directory with log files')
parser.add_argument('--no-plot', dest='no_plot', action='store_true', help='do not make histogram plots of the initial conditions; slightly improves speed')
parser.add_argument('--email', default='{email}', help='send Slurm email notifications to this address')
parser.add_argument('--nodes', type=int, default=20, help='how many nodes we want to request in Slurm submission script')
parser.add_argument('--tasks-per-node', dest='ntasks_per_node', type=int, default=40, help='how many MPI tasks per node we want to request in Slurm submission script')
parser.add_argument('--memory-per-task', dest='mem', default=4100, type=int, help='memory GIZMO will have available per task - inserted into .param file')
args = parser.parse_args()

fname_out = args.output_file
fnames_in = args.input_file
log_dir = args.log
seed = args.seed
copy = False if args.no_copy is True else True
load = args.load
no_plot = args.no_plot
email = args.email
nodes = args.nodes
ntasks_per_node = args.ntasks_per_node
memory_per_task = args.mem

if load is not None:
    if len(fnames_in) > 0:
        print('Do not specify input files if using --load.\n')
        parser.print_help()
        exit(1)
    elif seed is not None:
        print('Do not use both --seed and --load.\n')
        parser.print_help()
        exit(1)
    elif os.path.basename(os.path.realpath(fname_out)).startswith(os.path.basename(os.path.realpath(load))):
        print('Output file cannot be located in the directory we load from. This will cause errors if we try to load again in the future.')
        exit(1)
elif len(fnames_in) == 0:
    print('No input files specified! Stopping without creating any files.')
    exit(0)

# path to log directory
# for each galaxy the log directory will contain:
#   - a log file with position, velocity used, path given by get_log_file_path()
#   - a copy of .hdf5 data, path given by get_copied_data_file_path()
if log_dir is None:
    log_dir = fname_out + '.log'
elif os.path.basename(os.path.realpath(fname_out)).startswith(os.path.basename(os.path.realpath(log_dir))):
    print('Output file cannot be located in the log directory. This will cause errors if we try to load in the future.')
    exit(1)

# make the log directory if it doesn't already exist
try:
    os.mkdir(log_dir)
except FileExistsError:
    pass # if the directory already exists, do nothing

# if --load given, check that the directory from which we load exists
if load is not None:
    assert os.path.exists(load)

# initialize random number generator
if not load:
    if seed is None:
        seed = np.random.randint(2**63)
    random_gen = np.random.default_rng(seed)

# dtype for output
dtype = 'float64' if double_precision else 'float32'

# fields for all particles
common_fields = ['Coordinates', 'Masses', 'Velocities']
# fields for gas particles
gas_only_fields = ['Density', 'InternalEnergy'] if density_present else ['InternalEnergy']
# fields for star particles
star_only_fields = ['StellarFormationTime']
# fields for gas, star, disk, and bulge particles
gas_stars_only_fields = ['Metallicity'] if common.default_num_metals > 1 else []

gas_fields = (*common_fields, *gas_only_fields, *gas_stars_only_fields, )
dm_fields = (*common_fields, )
disk_fields = (*common_fields, )
bulge_fields = (*common_fields, )
star_fields = (*common_fields, *star_only_fields, *gas_stars_only_fields, )
bh_fields = (*common_fields, )

# a short, descriptive name for this set of initial conditions - e.g. this is used in file names
protocluster_name = os.path.basename(os.path.splitext(fname_out)[0])

# ----- functions -----
# additional functions defined in generate_protocluster_functions.py

# given a log directory and protocluster name, return path where the output plot with arrows is to be saved
def get_plot_path(log_dir, protocluster_name):
    return os.path.join(log_dir, protocluster_name + '.pdf')

# given a log directory and protocluster name, return path where the output plot without arrows is to be saved
def get_plot_path_no_arrows(log_dir, protocluster_name):
    return os.path.join(log_dir, protocluster_name + '_no_arrows.pdf')

# given a log directory and galaxy index, return path to that galaxy's log file
def get_log_file_path(log_dir, index):
    return os.path.join(log_dir, '{:0>3d}.pkl'.format(index))

# given a log directory and galaxy index, return path to that galaxy's .hdf5 file (inside the log directory)
def get_copied_data_file_path(log_dir, index):
    return os.path.join(log_dir, '{:0>3d}.hdf5'.format(index))

# given a log directory, return path where we will save .txt table of galaxy masses
def get_mass_table_path(log_dir):
    return os.path.join(log_dir, 'masses.txt')

# given a log directory and protocluster name, return path where we will save Slurm script
def get_slurm_script_path(log_dir, protocluster_name):
    return os.path.join(log_dir, protocluster_name + '.slurm_job.sh')

# given a log directory and protocluster name, return path where we will save .param GIZMO parameter file
def get_param_path(log_dir, protocluster_name):
    return os.path.join(log_dir, protocluster_name + '.param')


# ----- script part -----

# if --load was given, file names were not specified
# instead, we load them from the load directory (old log directory)
if load is not None:
    n_files = count_files_inside_with_suffix(load, '.hdf5')
    assert n_files > 0
    fnames_in = [get_copied_data_file_path(load, index) for index in range(n_files)]

# if --load was not given, initialize iterators that will give speeds and radial distances
else:
    # this iterator will yield distances from the centre of the protocluster for the galaxies
    # we can generate a random distance by calling next(it_r)
    it_r = sample_uniform(random_gen, b=r_maximum)

    # this iterator will yield speeds for the galaxies
    # we can generate a random speed by calling next(it_speed)
    it_speed = sample_half_gaussian_within_bounds(random_gen, v_std, v_maximum)

# these dictionaries will hold output - we initialize with empty lists that we'll append to
gas_lists = {field:[] for field in gas_fields}
dm_lists = {field:[] for field in dm_fields}
disk_lists = {field:[] for field in disk_fields}
bulge_lists = {field:[] for field in bulge_fields}
star_lists = {field:[] for field in star_fields}
bh_lists = {field:[] for field in bh_fields}

# these arrays will store the number of particles for each galaxy
Ngas = np.zeros(len(fnames_in), dtype='uint32')
Ndm = np.copy(Ngas)
Ndisk = np.copy(Ngas)
Nbulge = np.copy(Ngas)
Nstar = np.copy(Ngas)
Nbh = np.copy(Ngas)

# these arrays will store the total mass for each galaxy
Mgas = np.zeros(len(fnames_in), dtype='float64')
Mdm = np.copy(Mgas)
Mdisk = np.copy(Mgas)
Mbulge = np.copy(Mgas)
Mstar = np.copy(Mgas)
Mbh = np.copy(Mgas)

# these arrays will store the offset coordinates and bulk velocities of the galaxies
galaxy_coordinates = np.zeros((len(fnames_in), 3), dtype='float64')
galaxy_velocities = np.copy(galaxy_coordinates)

# these "empty" dictionaries will be used in place of bulge/star data if turned off
if not do_bulge:
    empty_data_bulge = {field:common.empty_array(field) for field in bulge_fields}
if not do_star:
    empty_data_star = {field:common.empty_array(field) for field in star_fields}

# iterate over all input files
for index, data_file in enumerate(fnames_in):
    print('Operating on galaxy {} at {}'.format(index, data_file))

    # convert data_file to absolute path (rather than a relative path)
    data_file = os.path.abspath(data_file)

    # copy/link input .hdf5 into the log directory
    # first check if we're loading from same place as our output - if so there's no need to copy
    if (load is None) or (os.path.realpath(log_dir) != os.path.realpath(load)):
        # case where we copy; at least one of:
        #   --no-copy was not given (so copy is True), or
        #   --no-copy and --load were both given and we're loading from a symbolic link (we'll copy it)
        if copy or ((load is not None) and (os.path.islink(data_file))):
            # the copy2 function preserves metadata
            shutil.copy2(data_file, get_copied_data_file_path(log_dir, index))
        # --no-copy given: use symbolic link with a relative path
        else: 
            link_name = get_copied_data_file_path(log_dir, index)
            os.symlink(os.path.relpath(data_file, start=os.path.basename(link_name)), link_name)            

    # load data from file
    with h5py.File(data_file, 'r') as f:
        this_gas = {field:common.grab_property(f, 0, field) for field in gas_fields}
        this_dm = {field:common.grab_property(f, 1, field) for field in dm_fields}
        this_disk = {field:common.grab_property(f, 2, field) for field in disk_fields}
        this_bulge = {field:common.grab_property(f, 3, field) for field in bulge_fields} if do_bulge else empty_data_bulge
        this_star = {field:common.grab_property(f, 4, field) for field in star_fields} if do_star else empty_data_star
        this_bh = {field:common.grab_property(f, 5, field) for field in bh_fields}

    # If we are using data for formed stars, we change their formation times
    # to negative times since the simulation will start at t=0.
    if do_star and len(this_star['StellarFormationTime']) > 0:
        this_star['StellarFormationTime'] = (this_star['StellarFormationTime']
                                             - 1.05 * np.max(this_star['StellarFormationTime']))

    Ngas[index] = len(this_gas['Coordinates'])
    Ndm[index] = len(this_dm['Coordinates'])
    Ndisk[index] = len(this_disk['Coordinates'])
    Nbulge[index] = len(this_bulge['Coordinates'])
    Nstar[index] = len(this_star['Coordinates'])
    Nbh[index] = len(this_bh['Coordinates'])

    Mgas[index] = np.sum(this_gas['Masses'])
    Mdm[index] = np.sum(this_dm['Masses'])
    Mdisk[index] = np.sum(this_disk['Masses'])
    Mbulge[index] = np.sum(this_bulge['Masses'])
    Mstar[index] = np.sum(this_star['Masses'])
    Mbh[index] = np.sum(this_bh['Masses'])

    if load is not None: # if --load was given, load the position and velocity for the galaxy
        with open(get_log_file_path(load, index), 'rb') as lf:
            old_data = pickle.load(lf)
        rand_vec = old_data['RandomPosition']
        rand_vel_vec = old_data['RandomVelocity']
        rotation = old_data['RandomRotationMatrix']
        assert index == old_data['Index']

    else: # if --load not given, make random position and velocity for galaxy
        # position offset of the entire galaxy 
        rand_vec = next(it_r) * rand_unit_vector(random_gen)
        # velocity offset for the entire galaxy
        rand_vel_vec = next(it_speed) * rand_unit_vector(random_gen)
        # rotation matrix that we will use to rotate the entire galaxy
        rotation = rand_rotation_matrix(random_gen)

    galaxy_coordinates[index,:] = rand_vec
    galaxy_velocities[index,:] = rand_vel_vec

    # save random position, velocity, rotation to the log, if we have not already loaded from it
    if (load is None) or (os.path.realpath(log_dir) != os.path.realpath(load)):
        log_data = {'Path':                 data_file,
                    'Index':                index,
                    'Seed':                 seed,
                    'RandomVelocity':       rand_vel_vec,
                    'RandomPosition':       rand_vec,
                    'RandomRotationMatrix': rotation}
        with open(get_log_file_path(log_dir, index), 'wb') as lf:
            pickle.dump(log_data, lf, protocol=pickle.HIGHEST_PROTOCOL)

    # adjust coordinates for new rotation and offset
    this_gas['Coordinates'] = set_new_coords(this_gas['Coordinates'], rotation, rand_vec)
    this_dm['Coordinates'] = set_new_coords(this_dm['Coordinates'], rotation, rand_vec)
    this_disk['Coordinates'] = set_new_coords(this_disk['Coordinates'], rotation, rand_vec)
    this_bulge['Coordinates'] = set_new_coords(this_bulge['Coordinates'], rotation, rand_vec)
    this_star['Coordinates'] = set_new_coords(this_star['Coordinates'], rotation, rand_vec)
    this_bh['Coordinates'] = set_new_coords(this_bh['Coordinates'], rotation, rand_vec)

    # adjust velocities for new rotation and offset
    this_gas['Velocities'] = set_new_coords(this_gas['Velocities'], rotation, rand_vel_vec)
    this_dm['Velocities'] = set_new_coords(this_dm['Velocities'], rotation, rand_vel_vec)
    this_disk['Velocities'] = set_new_coords(this_disk['Velocities'], rotation, rand_vel_vec)
    this_bulge['Velocities'] = set_new_coords(this_bulge['Velocities'], rotation, rand_vel_vec)
    this_star['Velocities'] = set_new_coords(this_star['Velocities'], rotation, rand_vel_vec)
    this_bh['Velocities'] = set_new_coords(this_bh['Velocities'], rotation, rand_vel_vec)

    # save data for this galaxy to the running lists
    for field in gas_fields:
        gas_lists[field].append(this_gas[field])
    for field in dm_fields:
        dm_lists[field].append(this_dm[field])
    for field in disk_fields:
        disk_lists[field].append(this_disk[field])
    for field in bulge_fields:
        bulge_lists[field].append(this_bulge[field])
    for field in star_fields:
        star_lists[field].append(this_star[field])
    for field in bh_fields:
        bh_lists[field].append(this_bh[field])

# combine data from all galaxies into arrays with all particles
gas = {field:np.concatenate(gas_lists[field]) for field in gas_fields}
dm = {field:np.concatenate(dm_lists[field]) for field in dm_fields}
disk = {field:np.concatenate(disk_lists[field]) for field in disk_fields}
bulge = {field:np.concatenate(bulge_lists[field]) for field in bulge_fields}
star = {field:np.concatenate(star_lists[field]) for field in star_fields}
bh = {field:np.concatenate(bh_lists[field]) for field in bh_fields}

# check that every field has an entry for each particle
assert all(len(gas[field]) == len(gas['Coordinates']) for field in gas_fields)
assert all(len(dm[field]) == len(dm['Coordinates']) for field in dm_fields)
assert all(len(disk[field]) == len(disk['Coordinates']) for field in disk_fields)
assert all(len(bulge[field]) == len(bulge['Coordinates']) for field in bulge_fields)
assert all(len(star[field]) == len(star['Coordinates']) for field in star_fields)
assert all(len(bh[field]) == len(bh['Coordinates']) for field in bh_fields)

# check that there are no bulge/star particles if we did not want them
if not do_bulge: assert len(bulge['Coordinates'])==0
if not do_star: assert len(star['Coordinates'])==0

# compute total number of particles
tot_Ngas = int(np.sum(Ngas))
tot_Ndm = int(np.sum(Ndm))
tot_Ndisk = int(np.sum(Ndisk))
tot_Nbulge = int(np.sum(Nbulge))
tot_Nstar = int(np.sum(Nstar))
tot_Nbh = int(np.sum(Nbh))
tot_part = int(tot_Ngas + tot_Ndm + tot_Ndisk + tot_Nbulge + tot_Nstar + tot_Nbh)
Npart = np.array([tot_Ngas, tot_Ndm, tot_Ndisk, tot_Nbulge, tot_Nstar, tot_Nbh], dtype='uint32')

# compute total mass (in code units)
tot_Mgas = np.sum(Mgas)
tot_Mdm = np.sum(Mdm)
tot_Mdisk = np.sum(Mdisk)
tot_Mbulge = np.sum(Mbulge)
tot_Mstar = np.sum(Mstar)
tot_Mbh = np.sum(Mbh)
tot_mass = tot_Mgas + tot_Mdm + tot_Mdisk + tot_Mbulge + tot_Mstar + tot_Mbh

# save table of galaxy masses to .txt file
# the `with` ignores division by zero warnings - e.g. this occurs if no bulge particle
with np.errstate(divide='ignore', invalid='ignore'):
    df_mass_table = DataFrame({
            'name':tuple((os.path.basename(os.path.splitext(fname)[0]) for fname in fnames_in)),
            'M':common.code_mass_to_Msun(Mgas+Mdm+Mdisk+Mbulge+Mstar+Mbh),
            'fgas':Mgas/(Mgas+Mdisk+Mbulge+Mstar),
            'Mgas':common.code_mass_to_Msun(Mgas),
            'Mdm':common.code_mass_to_Msun(Mdm),
            'Mdisk':common.code_mass_to_Msun(Mdisk),
            'Mbulge':common.code_mass_to_Msun(Mbulge),
            'Mstar':common.code_mass_to_Msun(Mstar),
            'Mbh':common.code_mass_to_Msun(Mbh),
            'Mstellar':common.code_mass_to_Msun(Mdisk+Mbulge+Mstar),
            'Mgas/particle':common.code_mass_to_Msun(Mgas/Ngas),
            'Mdm/particle':common.code_mass_to_Msun(Mdm/Ndm),
            'Mdisk/particle':common.code_mass_to_Msun(Mdisk/Ndisk),
            'Mbulge/particle':common.code_mass_to_Msun(Mbulge/Nbulge),
            'Mstar/particle':common.code_mass_to_Msun(Mstar/Nstar),
            'Mbh/particle':common.code_mass_to_Msun(Mbh/Nbh),
            'Mstellar/particle':common.code_mass_to_Msun((Mdisk+Mbulge+Mstar)/(Ndisk+Nbulge+Nstar))
        })
    df_mass_table.append({
            'name':'total',
            'M':common.code_mass_to_Msun(tot_mass),
            'fgas':tot_Mgas/(tot_Mgas+tot_Mdisk+tot_Mbulge+tot_Mstar),
            'Mgas':common.code_mass_to_Msun(tot_Mgas),
            'Mdm':common.code_mass_to_Msun(tot_Mdm),
            'Mdisk':common.code_mass_to_Msun(tot_Mdisk),
            'Mbulge':common.code_mass_to_Msun(tot_Mbulge),
            'Mstar':common.code_mass_to_Msun(tot_Mstar),
            'Mbh':common.code_mass_to_Msun(tot_Mbh),
            'Mstellar':common.code_mass_to_Msun(tot_Mdisk+tot_Mbulge+tot_Mstar),
            'Mgas/particle':common.code_mass_to_Msun(tot_Mgas/tot_Ngas),
            'Mdm/particle':common.code_mass_to_Msun(tot_Mdm/tot_Ndm),
            'Mdisk/particle':common.code_mass_to_Msun(tot_Mdisk/tot_Ndisk),
            'Mbulge/particle':common.code_mass_to_Msun(tot_Mbulge/tot_Nbulge),
            'Mstar/particle':common.code_mass_to_Msun(tot_Mstar/tot_Nstar),
            'Mbh/particle':common.code_mass_to_Msun(tot_Mbh/tot_Nbh),
            'Mstellar/particle':common.code_mass_to_Msun((tot_Mdisk+tot_Mbulge+tot_Mstar)
                                                         /(tot_Ndisk+tot_Nbulge+tot_Nstar))
        }, ignore_index=True)
df_mass_table.to_string(get_mass_table_path(log_dir), float_format='{:>15.7g}'.format)

# We need to remove the bulk velocity of the system. This is important so that
# the group of galaxies doesn't fly away rapidly.
#
# We start by mass weighting the velocities of all of the components, and then
# subtracting away the result of the sum.
print('Correcting for bulk velocities.')
mw_gv = np.sum(gas['Masses'].reshape((-1, 1)) * gas['Velocities'], axis=0)
mw_dmv = np.sum(dm['Masses'].reshape((-1, 1)) * dm['Velocities'], axis=0)
mw_dv = np.sum(disk['Masses'].reshape((-1, 1)) * disk['Velocities'], axis=0)
mw_bv = np.sum(bulge['Masses'].reshape((-1, 1)) * bulge['Velocities'], axis=0)
mw_sv = np.sum(star['Masses'].reshape((-1, 1)) * star['Velocities'], axis=0)
mw_bhv = np.sum(bh['Masses'].reshape((-1, 1)) * bh['Velocities'], axis=0)
bv = ((mw_gv + mw_dmv + mw_dv + mw_bv + mw_sv + mw_bhv) / tot_mass).reshape((1, 3))
print('Bulk velocity was: {}'.format(bv))
gas['Velocities'] -= bv
dm['Velocities'] -= bv
disk['Velocities'] -= bv
bulge['Velocities'] -= bv
star['Velocities'] -= bv
bh['Velocities'] -= bv
galaxy_velocities -= bv

# similarly, change coordinates so that the centre of mass is at (0, 0, 0)
print('Correcting for centre of mass.')
mw_gv = np.sum(gas['Masses'].reshape((-1, 1)) * gas['Coordinates'], axis=0)
mw_dmv = np.sum(dm['Masses'].reshape((-1, 1)) * dm['Coordinates'], axis=0)
mw_dv = np.sum(disk['Masses'].reshape((-1, 1)) * disk['Coordinates'], axis=0)
mw_bv = np.sum(bulge['Masses'].reshape((-1, 1)) * bulge['Coordinates'], axis=0)
mw_sv = np.sum(star['Masses'].reshape((-1, 1)) * star['Coordinates'], axis=0)
mw_bhv = np.sum(bh['Masses'].reshape((-1, 1)) * bh['Coordinates'], axis=0)
cm = ((mw_gv + mw_dmv + mw_dv + mw_bv + mw_sv + mw_bhv) / tot_mass).reshape((1, 3))
print('Centre of mass was: {}'.format(cm))
gas['Coordinates'] -= cm
dm['Coordinates'] -= cm
disk['Coordinates'] -= cm
bulge['Coordinates'] -= cm
star['Coordinates'] -= cm
bh['Coordinates'] -= cm
galaxy_coordinates -= cm

# plot the positions and velocities in 2D histograms with arrows
if not no_plot:
    print('Plotting.')
    fig, ax = plt.subplots(2, 3, figsize=(11, 8.5), dpi=600, constrained_layout=True)
    fig2, ax2 = plt.subplots(2, 3, figsize=(11, 8.5), dpi=600, constrained_layout=True)
    title = protocluster_name
    fig.suptitle(title, y=0.97, size='xx-large')
    fig2.suptitle(title, y=0.97, size='xx-large')
    # fig ====> will have arrows for velocities, galaxy labels      fig2 ====> won't have arrows, labels
    stellar_coords = np.concatenate(
        (disk['Coordinates'], bulge['Coordinates'], star['Coordinates']),
        axis=0)
    stellar_masses = np.concatenate(
        (disk['Masses'], bulge['Masses'], star['Masses']),
        axis=0)
    stellar_velocities = np.concatenate(
        (disk['Velocities'], bulge['Velocities'], star['Velocities']),
        axis=0)
    ax[0,0].set_xlabel('x [kpc]')
    ax[0,0].set_ylabel('z [kpc]')
    ax[1,0].set_xlabel('x [kpc]')
    ax[1,0].set_ylabel('z [kpc]')
    ax[0,1].set_xlabel('x [kpc]')
    ax[0,1].set_ylabel('y [kpc]')
    ax[1,1].set_xlabel('x [kpc]')
    ax[1,1].set_ylabel('y [kpc]')
    ax[0,2].set_xlabel('y [kpc]')
    ax[0,2].set_ylabel('z [kpc]')
    ax[1,2].set_xlabel('y [kpc]')
    ax[1,2].set_ylabel('z [kpc]')
    ax2[0,0].set_xlabel('x [kpc]')
    ax2[0,0].set_ylabel('z [kpc]')
    ax2[1,0].set_xlabel('x [kpc]')
    ax2[1,0].set_ylabel('z [kpc]')
    ax2[0,1].set_xlabel('x [kpc]')
    ax2[0,1].set_ylabel('y [kpc]')
    ax2[1,1].set_xlabel('x [kpc]')
    ax2[1,1].set_ylabel('y [kpc]')
    ax2[0,2].set_xlabel('y [kpc]')
    ax2[0,2].set_ylabel('z [kpc]')
    ax2[1,2].set_xlabel('y [kpc]')
    ax2[1,2].set_ylabel('z [kpc]')
    # plot 95 kpc circles for reference
    angle = np.linspace(0, 2*np.pi, num=1000)
    for r in [95]:
        for a in (*ax[1,:], *ax2[1,:]):
            a.plot(r*np.cos(angle), r*np.sin(angle), color='grey', linestyle=':')
            a.annotate(str(r) + ' kpc', np.array((-1,1))*np.sqrt(0.5)*(r+5), c='grey', rotation=45,
                       horizontalalignment='center', verticalalignment='center', size=9)
    print('Making histograms for upper subplots.')
    # plot zoomed in plots in the top row of subplots
    common.plot_2D_histogram(stellar_coords, cube_radius=75, axes=(*ax[0,:], *ax2[0,:]),
                             nbins=512, squish_along=[1,2,0,1,2,0], rasterized=True, nticks=7, 
                             mark_maximum=False, weights=stellar_masses, return_histogram=False)
    print('Making histograms for lower subplots.')
    # plot wider-field plots in the bottom row of subplots
    common.plot_2D_histogram(stellar_coords, cube_radius=110, axes=(*ax[1,:], *ax2[1,:]),
                             nbins=512, squish_along=[1,2,0,1,2,0], rasterized=True, nticks=6,
                             mark_maximum=False, weights=stellar_masses, return_histogram=False)
    # for each galaxy, plot an arrow on each axis to show the velocity
    arrow_length = 20./np.max(galaxy_velocities) # how long in kpc should velocity arrows be, per km/s
    # arrow properties
    ap = {'width':0.05, 'headwidth':0.4, 'headlength':0.4,
          'shrink':0, 'color':'chartreuse', 'alpha':0.7} 
    print('Plotting annotations.')
    for coord, vel, fname in zip(galaxy_coordinates, galaxy_velocities, fnames_in):
        # get a descriptive label for this galaxy (to include in plots) from the file name
        # are we loading and is fname a symbolic link?
        # if yes, get label from file name that link points to
        # if no, get label from fname itself
        label = (os.path.basename(os.path.splitext(os.readlink(fname))[0])
                 if (load is not None) and (os.path.islink(fname)) else
                 os.path.basename(os.path.splitext(fname)[0]))
        ax[0,0].annotate('', xy=coord[[0,2]]+arrow_length*vel[[0,2]],xytext=coord[[0,2]], arrowprops=ap)
        ax[1,0].annotate('', xy=coord[[0,2]]+arrow_length*vel[[0,2]],xytext=coord[[0,2]], arrowprops=ap)
        ax[0,1].annotate('', xy=coord[[0,1]]+arrow_length*vel[[0,1]],xytext=coord[[0,1]], arrowprops=ap)
        ax[1,1].annotate('', xy=coord[[0,1]]+arrow_length*vel[[0,1]],xytext=coord[[0,1]], arrowprops=ap)
        ax[0,2].annotate('', xy=coord[[1,2]]+arrow_length*vel[[1,2]],xytext=coord[[1,2]], arrowprops=ap)
        ax[1,2].annotate('', xy=coord[[1,2]]+arrow_length*vel[[1,2]],xytext=coord[[1,2]], arrowprops=ap)
        ax[0,0].annotate(label, xy=coord[[0,2]], fontsize=8)
        ax[1,0].annotate(label, xy=coord[[0,2]], fontsize=8)
        ax[0,1].annotate(label, xy=coord[[0,1]], fontsize=8)
        ax[1,1].annotate(label, xy=coord[[0,1]], fontsize=8)
        ax[0,2].annotate(label, xy=coord[[1,2]], fontsize=8)
        ax[1,2].annotate(label, xy=coord[[1,2]], fontsize=8)
    # save figure to file
    print('Saving plots.')
    fig.savefig(get_plot_path(log_dir, protocluster_name))
    fig2.savefig(get_plot_path_no_arrows(log_dir, protocluster_name))
    # delete the figure to free up memory
    print('Clearing figures.')
    fig.clf()
    fig2.clf()
    plt.close(fig)
    plt.close(fig2)

# generate HDF5 file
print('Generating HDF5 file.')

# generate new particle IDs, starting at 1
new_ids = np.arange(1, tot_part + 1, 1, dtype='uint32')

with h5py.File(fname_out, 'w') as fp:

    # set header properties
    header = fp.create_group('Header')
    header.attrs['MassTable'] = np.zeros(6);
    header.attrs['Time'] = 0.0;  # initial time
    # header.attrs['Redshift'] = 0.0; # initial redshift
    # header.attrs['BoxSize'] = 1.0e4; # box size
    header.attrs['NumFilesPerSnapshot'] = 1; # number of files 
    # header.attrs['Omega0'] = 0.; # z=0 Omega_matter
    # header.attrs['OmegaLambda'] = 0.; # z=0 Omega_Lambda
    # header.attrs['HubbleParam'] = 1.0; # z=0 hubble parameter (small 'h'=H/100 km/s/Mpc)
    # header.attrs['Flag_Sfr'] = 1; # flag indicating whether star formation is on or off
    # header.attrs['Flag_Cooling'] = 1; # flag indicating whether cooling is on or off
    # header.attrs['Flag_StellarAge'] = 1; # flag indicating whether stellar ages are to be saved
    # header.attrs['Flag_Metals'] = 11; # flag indicating whether metallicity are to be saved
    # header.attrs['Flag_Feedback'] = 1; # flag indicating whether some parts of springel-hernquist model are active
    header.attrs['Flag_DoublePrecision'] = double_precision; # flag indicating whether ICs are in single/double precision
    # header.attrs['Flag_IC_Info'] = 0; # flag indicating extra options for ICs
    header.attrs['NumPart_ThisFile'] = Npart;
    header.attrs['NumPart_Total'] = Npart;
    header.attrs['NumPart_Total_HighWord'] = 0 * Npart;

    # running index for the new ParticleIDs
    running_idx = 0

    # gas
    if len(gas['Coordinates']) > 0:
        print('Adding gas particles.')
        p = fp.create_group('PartType0')
        p.create_dataset('ParticleIDs', data=new_ids[running_idx:running_idx + tot_Ngas], dtype='uint32')
        for field in gas_fields:
            p.create_dataset(field, data=gas[field], dtype=dtype)
        running_idx += tot_Ngas
    else:
        print('No gas particles.')

    # dark matter
    if len(dm['Coordinates']) > 0:
        print('Adding dark matter particles.')
        p = fp.create_group('PartType1')
        p.create_dataset('ParticleIDs', data=new_ids[running_idx:running_idx + tot_Ndm], dtype='uint32')
        for field in dm_fields:
            p.create_dataset(field, data=dm[field], dtype=dtype)
        running_idx += tot_Ndm
    else:
        print('No dark matter particles.')

    # disk
    if len(disk['Coordinates']) > 0:
        print('Adding disk particles.')
        p = fp.create_group('PartType2')
        p.create_dataset('ParticleIDs', data=new_ids[running_idx:running_idx + tot_Ndisk], dtype='uint32')
        for field in disk_fields:
            p.create_dataset(field, data=disk[field], dtype=dtype)
        running_idx += tot_Ndisk
    else:
        print('No disk particles.')

    # bulge
    if len(bulge['Coordinates']) > 0:
        print('Adding bulge particles.')
        p = fp.create_group('PartType3')
        p.create_dataset('ParticleIDs', data=new_ids[running_idx:running_idx + tot_Nbulge], dtype='uint32')
        for field in bulge_fields:
            p.create_dataset(field, data=bulge[field], dtype=dtype)
        running_idx += tot_Nbulge
    else:
        print('No bulge particles.')

    # stars formed
    if len(star['Coordinates']) > 0:
        print('Adding star particles.')
        p = fp.create_group('PartType4')
        p.create_dataset('ParticleIDs', data=new_ids[running_idx:running_idx + tot_Nstar], dtype='uint32')
        for field in star_fields:
            p.create_dataset(field, data=star[field], dtype=dtype)
        running_idx += tot_Nstar
    else:
        print('No star particles.')

    # black holes
    if len(bh['Coordinates']) > 0:
        print('Adding black hole particles.')
        p = fp.create_group('PartType5')
        p.create_dataset('ParticleIDs', data=new_ids[running_idx:running_idx + tot_Nbh], dtype='uint32')
        for field in bh_fields:
            p.create_dataset(field, data=bh[field], dtype=dtype)
        running_idx += tot_Nbh
    else:
        print('No black hole particles.')

    # check that we've gone through all of the particles
    assert running_idx == tot_part

print('Done!')


# generate GIZMO parameter file and Slurm job script for these initial conditions

# load the templates
with open(path_to_slurm_script_template, 'r') as f:
    sh_template = f.read();
with open(path_to_gizmo_param_template, 'r') as f:
    param_template = f.read()

# formatter maps fields to the values we use (with a __missing__ method for handling fields not set) 

formatter = {
        'ic_file_name':os.path.basename(os.path.splitext(fname_out)[0]), # we want no extension for this
        'param_file_name':os.path.basename(get_param_path(log_dir, protocluster_name)),
        'job_name':protocluster_name,
        'resubmit_slurm_script':os.path.basename(get_slurm_script_path(log_dir, protocluster_name)),
        'timelim_string':'1-00:00',
        'timelim_in_s':'86400',
        'nodes':nodes,
        'ntasks_per_node':ntasks_per_node,
        'memory_per_node':'0', # setting this to 0 reserves all memory on the node
        'memory_per_task':'{:d}'.format(memory_per_task),
        'resubmit_on':'1',
        'sim_time_max':'2.5',
        'black_hole_res':'{:.6g}'.format(max(tot_Mbh/tot_Nbh, tot_Mdm/tot_Ndm)), # should be at least as large as halo resolution
        'gas_fraction':'{:.6g}'.format(tot_Mgas/(tot_Mgas+tot_Mdisk+tot_Mbulge+tot_Mstar)),
        'stellar_mass':'{:.6g}'.format(common.code_mass_to_Msun(tot_Mdisk+tot_Mbulge+tot_Mstar)),
        'halo_mass':'{:.6g}'.format(common.code_mass_to_Msun(tot_Mdm)),
        'email':email,
        'gizmo_commit':'d768f3d'
    }
with open(get_slurm_script_path(log_dir, protocluster_name), 'w') as f:
    f.write(sh_template.format_map(formatter))
with open(get_param_path(log_dir, protocluster_name), 'w') as f:
    f.write(param_template.format_map(formatter))


