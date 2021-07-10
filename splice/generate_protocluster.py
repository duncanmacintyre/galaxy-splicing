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

# ----- constants -----

double_precision = True  # use True if GIZMO expects double precision ICs, False if single precision
do_bulge = False         # whether to include bulge particles in output (we always include disk particles)
do_star = False          # whether to include formed star particles in output (we always include disk particles)
density_present = False  # whether to include the gas field Density in output
num_metals = 11          # how many metals we track - set to 0 to ignore Metallicity fields
r_std = 59.449           # standard deviation for galaxy distance from centre, in kpc - defines normal distribution
r_maximum = 138.88       # maximum allowed galaxy distance from centre, in kpc
v_std = 571.83           # standard deviation for galaxy speeds, in km/s - defines normal distribution
v_maximum = 1681.9       # maximum allowed galaxy speed, in km/s

# these values for r_std, r_maximum, v_std, and v_maximum are
# derived in the file ../generate/position_velocity_distribution.ipynb

# ----- setup -----

parser = argparse.ArgumentParser(description='Splice together individual galaxies to form a cluster.\n Combines .hdf5 GIZMO snapshot files to generate one .hdf5 initial conditions file.',)
parser.add_argument('output_file', type=str, help='path to output HDF5 file to generate')
parser.add_argument('input_file', nargs='*', type=str, help='paths to HDF5 files to splice together')
parser.add_argument('--seed', type=int, default=None, help='integer seed for generating random variables; omit for random seed')
parser.add_argument('--no-copy', action='store_true', dest='no_copy', help="don't copy the input files into the log folder, use symbolic links instead")
parser.add_argument('--load', metavar='DIR', help='load a previously stored configuration; DIR is the directory with log files; WARNING: --load not fully tested, could be buggy')
parser.add_argument('--mass-table', action='store_true', dest='mass_table', help='for input files, use MassTable from headers; otherwise we get mass for each particle using Masses fields')
args = parser.parse_args()

fname_out = args.output_file
fnames_in = args.input_file
seed = args.seed
no_copy = args.no_copy
load = args.load
mass_table = args.mass_table

if load is not None:
    if len(fnames_in) > 0:
        print('Do not specify input files if using --load.\n')
        parser.print_help()
        exit(1)
    elif seed is not None:
        print('Do not use both --seed and --load.\n')
        parser.print_help()
        exit(1)
elif len(fnames_in) == 0:
    print('No input files specified! Stopping without creating any files.')
    exit(0)

# path to log directory
# for each galaxy the log directory will contain:
#   - a log file with position, velocity used, path given by get_log_file_path()
#   - a copy of .hdf5 data, path given by get_copied_data_file_path()
log_dir = fname_out + '.log'
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
# fields for gas, star, disk, and bulge particles
gas_stars_only_fields = ['Metallicity'] if num_metals > 1 else []

# how many items per field
field_size = {
    'Coordinates':     3,
    'Masses':          1,
    'Velocities':      3,
    'Density':         1,
    'InternalEnergy':  1,
    'Metallicity':     num_metals
}

gas_fields = (*common_fields, *gas_only_fields, *gas_stars_only_fields, )
dm_fields = (*common_fields, )
disk_fields = (*common_fields, )
bulge_fields = (*common_fields, )
star_fields = (*common_fields, *gas_stars_only_fields, )
bh_fields = (*common_fields, )


# ----- functions -----

def rand_rotation_matrix(deflection = 1.0, randnums = None):
    """
    Creates a random rotation matrix.
    
    deflection: the magnitude of the rotation. For 0, no rotation; for 1, 
    competely random rotation. Small deflection => small perturbation.
    randnums: 3 random numbers in the range [0, 1]. If `None`, 
    they will be auto-generated.
    """
    # from http://www.realtimerendering.com/resources/GraphicsGems/
    #gemsiii/rand_rotation.c
    
    if randnums is None:
        randnums = random_gen.uniform(size = (3,))
        
    theta, phi, z = randnums
    
    theta = theta * 2.0 * deflection * np.pi  # Rotation about the pole (Z).
    phi = phi * 2.0 * np.pi  # For direction of pole deflection.
    z = z * 2.0 * deflection  # For magnitude of pole deflection.
    
    # Compute a vector V used for distributing points over the sphere
    # via the reflection I - V Transpose(V).  This formulation of V
    # will guarantee that if x[1] and x[2] are uniformly distributed,
    # the reflected points will be uniform on the sphere.  Note that V
    # has length sqrt(2) to eliminate the 2 in the Householder matrix.
    
    r = np.sqrt(z)
    V = (
            np.sin(phi) * r,
            np.cos(phi) * r,
            np.sqrt(2.0 - z)
        )
    
    st = np.sin(theta)
    ct = np.cos(theta)
    
    R = np.array(((ct, st, 0), (-st, ct, 0), (0, 0, 1)))
    
    # Construct the rotation matrix  ( V Transpose(V) - I ) R.
    
    M = (np.outer(V, V) - np.eye(3)).dot(R)
    return M


# return a length-3 unit vector in a random direction (all directions equally likely)
def rand_unit_vector():
    phi = random_gen.uniform(0, 2.*np.pi)
    costheta = random_gen.uniform(-1, 1)
    sintheta = np.sin(np.arccos(costheta))
    return np.array((sintheta * np.cos(phi), sintheta * np.sin(phi), costheta))


# GENERATOR: sample half-normal distribution with mode at 0 and yield only values <= maximum
def sample_half_gaussian_within_bounds(std, maximum):
    while True:
        mag = np.abs(random_gen.normal(0, std, 1))
        if mag <= maximum:
            yield mag


# given Nx3 array, apply rotation matrix to each row then translate by vector offset
def set_new_coords(coords_vec, rotation, translation):
    for i, vec in enumerate(coords_vec):
        coords_vec[i, :] = np.dot(rotation, vec)
    return coords_vec + np.asarray(translation).reshape((1, 3))


# given h5py.File f, return values for a part type and field or an empty array if not present
def grab_property(f, part_type, field):
    try:
        prop = np.asarray(f['/PartType%d/%s' % (part_type, field)])
    except KeyError:
        prop = empty_array(field)
        #print('KeyError: PartType%d/%s' % (part_type, field))
    return prop


# return empty numpy array of suitable size for given field
def empty_array(field):
    n = field_size[field] # how many entries per particle
    return np.empty((0,) if n==1 else (0, n))

# return how many files inside a given directory have names ending with suffix
def count_files_inside_with_suffix(path, suffix):
    return sum(1 for _ in filter(lambda s: s.endswith(suffix), os.listdir(path)))


# given a log directory and galaxy index, return path to that galaxy's log file
def get_log_file_path(log_dir, index):
    return os.path.join(log_dir, '{:0>3d}.pkl'.format(index))


# given a log directory and galaxy index, return path to that galaxy's .hdf5 file (inside the log directory)
def get_copied_data_file_path(log_dir, index):
    return os.path.join(log_dir, '{:0>3d}.hdf5'.format(index))


# ----- script part -----

# this iterator will yield distances from the centre of the protocluster for the galaxies
# we can generate a random distance by calling next(it_r)
it_r = sample_half_gaussian_within_bounds(r_std, r_maximum)

# this iterator will yield speeds for the galaxies
# we can generate a random speed by calling next(it_speed)
it_speed = sample_half_gaussian_within_bounds(v_std, v_maximum)

# if --load was given, file names were not specified
# instead, we load them from the load directory (old log directory)
if load is not None:
    n_files = count_files_inside_with_suffix(load, '.hdf5')
    assert n_files > 0
    fnames_in = [get_copied_data_file_path(load, index) for index in range(n_files)]

# these dictionaries will hold output - we initialize with empty lists that we'll append to
gas_lists = {field:[] for field in gas_fields}
dm_lists = {field:[] for field in dm_fields}
disk_lists = {field:[] for field in disk_fields}
bulge_lists = {field:[] for field in bulge_fields}
star_lists = {field:[] for field in star_fields}
bh_lists = {field:[] for field in bh_fields}

# these arrays will store the number of particles
Ngas = np.zeros(len(fnames_in), dtype = 'uint32')
Ndm = np.copy(Ngas)
Ndisk = np.copy(Ngas)
Nbulge = np.copy(Ngas)
Nstar = np.copy(Ngas)
Nbh = np.copy(Ngas)

# these "empty" dictionaries will be used in place of bulge/star data if turned off
if not do_bulge:
    empty_data_bulge = {field:empty_array(field) for field in bulge_fields}
if not do_star:
    empty_data_star = {field:empty_array(field) for field in star_fields}

# iterate over all input files
for index, data_file in enumerate(fnames_in):
    print('Operating on galaxy {} at {}'.format(index, data_file))

    # convert data_file to absolute path
    data_file = os.path.abspath(data_file)

    # copy/link input .hdf5 into the log directory
    if no_copy: # --no-copy given: use symbolic link
        os.symlink(data_file, get_copied_data_file_path(log_dir, index))
    else: # --no-copy wasn't given: copy the file
        # the copy2 function preserves metadata
        shutil.copy2(data_file, get_copied_data_file_path(log_dir, index))

    # load data from file
    with h5py.File(data_file, 'r') as f:
        this_gas = {field:grab_property(f, 0, field) for field in gas_fields}
        this_dm = {field:grab_property(f, 1, field) for field in dm_fields}
        this_disk = {field:grab_property(f, 2, field) for field in disk_fields}
        this_bulge = {field:grab_property(f, 3, field) for field in bulge_fields} if do_bulge else empty_data_bulge
        this_star = {field:grab_property(f, 4, field) for field in star_fields} if do_star else empty_data_star
        this_bh = {field:grab_property(f, 5, field) for field in bh_fields}

    Ngas[index] = len(this_gas['Coordinates'])
    Ndm[index] = len(this_dm['Coordinates'])
    Ndisk[index] = len(this_disk['Coordinates'])
    Nbulge[index] = len(this_bulge['Coordinates'])
    Nstar[index] = len(this_star['Coordinates'])
    Nbh[index] = len(this_bh['Coordinates'])

    # We don't use /PartTypeX/Masses if mass_table is True. In this case
    # we have to set them manually using the MassTable header data.
    if mass_table:
        mass_table_values = f['/Header'].attrs['MassTable']
        
        this_gas['Masses'] = np.ones(Ngas[index])
        this_dm['Masses'] = np.ones(Ndm[index])
        this_disk['Masses'] = np.ones(Ndisk[index])
        this_bulge['Masses'] = np.ones(Nbulge[index])
        this_star['Masses'] = np.ones(Nstar[index])
        this_bh['Masses'] = np.ones(Nbh[index])

        this_gas['Masses'][:] = mass_table_values[0]
        this_dm['Masses'][:] = mass_table_values[1]
        this_disk['Masses'][:] = mass_table_values[2]
        this_bulge['Masses'][:] = mass_table_values[3]
        this_star['Masses'][:] = mass_table_values[4]
        this_bh['Masses'][:] = mass_table_values[5]


    if load is not None: # if --load was given, load the position and velocity for the galaxy
        with open(get_log_file_path(load, index), 'rb') as lf:
            old_data = pickle.load(lf)
        rand_vec = old_data['RandomPosition']
        rand_vel_vec = old_data['RandomVelocity']
        rotation = old_data['RandomRotationMatrix']
        assert index == old_data['Index']

    else: # if --load not given, make random position and velocity for galaxy
        # position offset of the entire galaxy 
        rand_vec = next(it_r) * rand_unit_vector()
        # velocity offset for the entire galaxy
        rand_vel_vec = next(it_speed) * rand_unit_vector()
        # rotation matrix that we will use to rotate the entire galaxy
        rotation = rand_rotation_matrix()

    # save random position, velocity, rotation to the log
    log_data = {'Path':                 data_file,
                'Index':                index,
                'Seed':                 seed,
                'RandomVelocity':       rand_vel_vec,
                'RandomPosition':       rand_vec,
                'RandomRotationMatrix': rotation}
    with open(get_log_file_path(log_dir, index), 'wb') as lf:
        pickle.dump(log_data, lf, protocol = pickle.HIGHEST_PROTOCOL)

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

# compute total mass
tot_gas_mass = np.sum(gas['Masses'])
tot_dm_mass = np.sum(dm['Masses'])
tot_disk_mass = np.sum(disk['Masses'])
tot_bulge_mass = np.sum(bulge['Masses'])
tot_star_mass = np.sum(star['Masses'])
tot_bh_mass = np.sum(bh['Masses'])
tot_mass = tot_gas_mass + tot_dm_mass + tot_disk_mass + tot_bulge_mass + tot_star_mass + tot_bh_mass

# compute total number of particles
tot_Ngas = int(np.sum(Ngas))
tot_Ndm = int(np.sum(Ndm))
tot_Ndisk = int(np.sum(Ndisk))
tot_Nbulge = int(np.sum(Nbulge))
tot_Nstar = int(np.sum(Nstar))
tot_Nbh = int(np.sum(Nbh))
tot_part = int(tot_Ngas + tot_Ndm + tot_Ndisk + tot_Nbulge + tot_Nstar + tot_Nbh)
Npart = np.array([tot_Ngas, tot_Ndm, tot_Ndisk, tot_Nbulge, tot_Nstar, tot_Nbh], dtype='uint32')

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
print('Bulk velocity: {}'.format(bv))
gas['Velocities'] -= bv
dm['Velocities'] -= bv
disk['Velocities'] -= bv
bulge['Velocities'] -= bv
star['Velocities'] -= bv
bh['Velocities'] -= bv

# generate new particle IDs, starting at 1
new_ids = np.arange(1, tot_part + 1, 1, dtype='uint32')

# generate HDF5 file

print('Generating HDF5 file.')
with h5py.File(fname_out, 'w') as fp:

    # !!! verify which properties we need - see manual
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


