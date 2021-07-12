#!/bin/python
# DR: original version of code for stellar mass (for 2020 Rennehan et al. paper)
# DM: vastly rewrote to make it easy to change radii, added more output files,
#     refactored for speed, memory, readability, and versatility (June 2021)
#
# use --help for usage
# warning: currently there will be bugs if --mmg is used
#

import h5py
import numpy as np
from argparse import ArgumentParser
import os.path
import sys
import multiprocessing

from common import locate_peak_density_3D

# ---------- argument parsing

parser = ArgumentParser(description='Extract masses from GIZMO snapshot files, optionally computing mass contained in radial bins')
parser.add_argument('data_dir', help='path to directory containing snapshot files')
g1 = parser.add_mutually_exclusive_group(required=True)
g1.add_argument('-R', nargs='+', type=float, help='compute mass within bins - will use radii R [R [R ...]]')
g1.add_argument('-r', nargs=2, metavar=('DR', 'MAX_R'), type=float, help='compute mass within bins - will use radii DR, 2DR, 3DR, ... up to but not including MAX_R')
g1.add_argument('--total-only', action='store_true', dest='total_only', help='compute total masses only, do not bin by radii')
g2 = parser.add_mutually_exclusive_group(required=True)
g2.add_argument('-N', nargs='+', type=int, help='the snapshot file numbers to use')
g2.add_argument('-n', type=int, metavar='MAX_SNAP', help='use snapshot numbers 0, 1, 2, 3, ... up to but not including MAX_SNAP; e.g. use 100 for 0 through 99')
g3 = parser.add_mutually_exclusive_group(required=False)
g3.add_argument('-T', nargs='+', type=float, help='the times for each snapshot number')
g3.add_argument('-t', metavar='DT', type=float, help='snapshots have time N*DT where N is snapshot number')
g4 = parser.add_mutually_exclusive_group(required=False)
g4.add_argument('--print', action='store_true', help='print the results instead of saving them to files')
g4.add_argument('--prefix', help='prepend PREFIX to all output file names')
parser.add_argument('--no-weight', action='store_true', dest='no_weight', help='when binning, set origin at max star particle density rather than max star mass density; improves speed by about an order of magnitude')
parser.add_argument('--where', action='store_true', help='use numpy.where and matrix sum instead of numpy.fromiter and iterative sum; may or may not be faster')
parser.add_argument('--mmg', action='store_true', help='use most massive galaxy for centre rather than peak brightness; not currently working')
parser.add_argument('-d', default='  ', help='separator between fields in the output files; default two spaces')
parser.add_argument('-p', default=int(6), type=int, help='how many significant figures to use in the output file; default 6')
parser.add_argument('-w', default=int(11), type=int, help='all fields padded to be at least this width in the output files; default 11')
args = parser.parse_args()
    
# data_dir - the path to the folder containing snapshot files
data_dir = args.data_dir

# R - the radii we will use
if args.total_only: # --total-only was given
    R = np.array([]) # we are not binning by radii
elif args.R is not None: # -R was given
    R = np.array(args.R)
    assert(len(R) > 0)
elif args.r is not None: # -r was given
    R = np.arange(args.r[0], args.r[1], args.r[0])
    assert(len(R) > 0)
else:
    assert(False) # this should never be reached
assert(all(R > 0))

# snaps - the numbers of the snapshot files to process
if args.N is not None: # -N was given
    snaps = np.array(args.N)
else: # -n was given
    snaps = np.arange(0, args.n)
assert(len(snaps) > 0)
assert(all(snaps >= 0))

# t - the time associated with snapshots
if args.T is not None:
    t = np.array(args.T)
    assert(all(t >= 0))
    assert(len(t) == len(snaps))
elif args.t is not None:
    assert(args.t > 0)
    t = (args.t * snaps)
else:
    t = None

# output formatting
delimiter = args.d # separator between fields in the output file
precision = args.p # significant figures to use in the output file
assert(precision > 0)
fieldwidth = args.w # all fields padded to be at least this width in the output file
assert(fieldwidth > 0)


# ---------- functions

# --- formatter
# convert float to string; for outputs
formatter = ('{:>' + str(fieldwidth) + '.' + str(precision) + 'g}').format

# --- snap_to_fname
# given a snapshot number, return its file name (without full path)
def snap_to_fname(snap):
    return 'snapshot_{}.hdf5'.format(str(snap).zfill(3))

# --- get_mass_in_cylinder
# given particle masses and coords, compute stellar mass contained within radii in R
#   column_to_remove=None         sphere
#   column_to_remove=0            cylinder about x-axis
#   column_to_remove=1            cylinder about y-axis
#   column_to_remove=2            cylinder about z-axis
#
#   masses, coords should by numpy arrays
#
if not args.where: # --where was not specified
    def get_mass_in_cylinder(masses, coords, column_to_remove):
        if column_to_remove is None:
            radii = np.sqrt(coords[:,0]**2 + coords[:,1]**2 + coords[:,2]**2)
        elif column_to_remove == 0:
            radii = np.sqrt(coords[:,1]**2 + coords[:,2]**2)
        elif column_to_remove == 1:
            radii = np.sqrt(coords[:,0]**2 + coords[:,2]**2)
        elif column_to_remove == 2:
            radii = np.sqrt(coords[:,0]**2 + coords[:,1]**2)
        else:
            assert(False) # this line should never be reached

        return np.fromiter((masses[radii < i].sum() for i in R), count=len(R), dtype='float')
else: # --where was specified
    def get_mass_in_cylinder(masses, coords, column_to_remove):
        if column_to_remove is None:
            radii = np.sqrt(coords[:,0]**2 + coords[:,1]**2 + coords[:,2]**2)
        elif column_to_remove == 0:
            radii = np.sqrt(coords[:,1]**2 + coords[:,2]**2)
        elif column_to_remove == 1:
            radii = np.sqrt(coords[:,0]**2 + coords[:,2]**2)
        elif column_to_remove == 2:
            radii = np.sqrt(coords[:,0]**2 + coords[:,1]**2)
        else:
            assert(False) # this line should never be reached

        return np.where(radii.reshape(1,-1) < R.reshape(-1,1), masses.reshape(1,-1), 0).sum(axis=1)


# --- recenter
# given masses and coordinates, return coordinates shifted to center on peak brightness or most massive galaxy
if args.mmg: # --mmg was specified
    particle_ids_center = np.loadtxt(os.path.join(data_dir, 'mmg_particle_ids.txt'))
    def recenter(masses, coords, disk_ids):   # this function won't currently be working as written
        if data_dir == 'SPT2349_1e6_gf0.9':
            mmg_idx = np.array([not i for i in np.in1d(disk_ids, particle_ids_center)])
        else:
            mmg_idx = np.in1d(disk_ids, particle_ids_center)
        masses = masses[mmg_idx]
        coords = coords[mmg_idx]

        mmg_mass_sum = np.sum(mmg_disk_masses)
        mmg_com_x = np.sum(mmg_disk_masses * mmg_disk_coords[:, 0]) / mmg_mass_sum
        mmg_com_y = np.sum(mmg_disk_masses * mmg_disk_coords[:, 1]) / mmg_mass_sum
        mmg_com_z = np.sum(mmg_disk_masses * mmg_disk_coords[:, 2]) / mmg_mass_sum
        
        center_offset = np.array([mmg_com_x, mmg_com_y, mmg_com_z])
        
        # center on the peak brightness or most massive galaxy
        return coords - center_offset

# --- process_snapshot
# given one snapshot, return a tuple of fours strings giving the mass-binning output
#
#   snap is the snapshot number
#
# returned strings are for total stars, stars formed in sim, gas, and halo respectively
# each output string has 2+5*len(R) items
#   first item: total mass
#   second item: how many particles
#   next len(R) items: mass within spheres of radii from R
#   next len(R) items: mass within cylinders about x-axis of radii from R
#   next len(R) items: mass within cylinders about y-axis of radii from R
#   next len(R) items: mass within cylinders about z-axis of radii from R
#   next len(R) items: mass within cylinders, mean for the three orientations
# if time is not None, will have an additional item for time first (3+5*len(R) in total)
# each item is fieldwidth long (or longer, if precision too high); items are separated by delimiter
#
if not args.mmg: # --mmg was not specified
    def process_snapshot(snap, time=None):

        data_file = os.path.join(data_dir, 'snapshot_{}.hdf5'.format(str(snap).zfill(3)))
        print('Operating on {}'.format(data_file))
        sys.stdout.flush()

        with h5py.File(data_file, 'r') as f:
            disk_masses = np.array(f['/PartType2/Masses']) * 1e10
            disk_coords = np.array(f['/PartType2/Coordinates'])
            try: # it may be that no stars have formed yet, in which case we get an error
                formed_stellar_masses = np.array(f['/PartType4/Masses']) * 1e10
                formed_stellar_coords = np.array(f['/PartType4/Coordinates'])
            except:
                formed_stellar_masses = np.empty((0,), dtype=disk_masses.dtype)
                formed_stellar_coords = np.empty((0, 3), dtype=disk_coords.dtype)
            gas_masses = np.array(f['/PartType0/Masses']) * 1e10
            gas_coords = np.array(f['/PartType0/Coordinates'])
            halo_masses = np.array(f['/PartType1/Masses']) * 1e10
            halo_coords = np.array(f['/PartType1/Coordinates'])
            
        # all stars, both initial and formed
        star_masses = np.concatenate((disk_masses, formed_stellar_masses), axis=0)
        star_coords = np.concatenate((disk_coords, formed_stellar_coords), axis=0)

        # find peak brightness
        if args.no_weight:
            centre = locate_peak_density_3D(
                star_coords, cube_radius=75, nbins=512).reshape((1, 3))
        else:
            centre = locate_peak_density_3D(
                star_coords, cube_radius=75, nbins=512, weights=star_masses).reshape((1, 3))
        return tuple(process_snapshot_helper(m, c - centre, time) for m, c in (
                (star_masses, star_coords),
                (formed_stellar_masses, formed_stellar_coords),
                (gas_masses, gas_coords),
                (halo_masses, halo_coords),
            ))

else: # --mmg was specified
    def process_snapshot(snap, time=None):

        data_file = os.path.join(data_dir, 'snapshot_{}.hdf5'.format(str(snap).zfill(3)))
        print('Operating on {}'.format(data_file))
        sys.stdout.flush()

        with h5py.File(data_file, 'r') as f:
            disk_ids = np.array(f['/PartType2/ParticleIDs'])
            disk_masses = np.array(f['/PartType2/Masses']) * 1e10
            disk_coords = np.array(f['/PartType2/Coordinates'])
            try: # it may be that no stars have formed yet, in which case we get an error
                formed_stellar_masses = np.array(f['/PartType4/Masses']) * 1e10
                formed_stellar_coords = np.array(f['/PartType4/Coordinates'])
            except:
                formed_stellar_masses = np.empty((0, 1), dtype=disk_masses.dtype)
                formed_stellar_coords = np.empty((0, 3), dtype=disk_coords.dtype)
            gas_masses = np.array(f['/PartType0/Masses']) * 1e10
            gas_coords = np.array(f['/PartType0/Coordinates'])
            halo_masses = np.array(f['/PartType1/Masses']) * 1e10
            halo_coords = np.array(f['/PartType1/Coordinates'])
            
            # all stars, both initial and formed
            star_masses = np.concatenate((disk_masses, formed_stellar_masses), axis=0)
            star_coords = np.concatenate((disk_coords, formed_stellar_coords), axis=0)

        return tuple(process_snapshot_helper(m, recenter(m, c, disk_ids), time) for m, c in (
                (star_masses, star_coords),
                (formed_stellar_masses, formed_stellar_coords),
                (gas_masses, gas_coords),
                (halo_masses, halo_coords),
            ))

# given masses and coords, bin & sum masses and return the results string
def process_snapshot_helper(masses, coords, time):
    s = get_mass_in_cylinder(masses, coords, None)
    x = get_mass_in_cylinder(masses, coords, 0)
    y = get_mass_in_cylinder(masses, coords, 1)
    z = get_mass_in_cylinder(masses, coords, 2)
    mean = (x + y + z) / 3.
    if time is not None:
        return delimiter.join(formatter(x) for x in (time, np.sum(masses), len(masses), *s, *x, *y, *z, *z, *mean))
    else:
        return delimiter.join(formatter(x) for x in (np.sum(masses), len(masses), *s, *x, *y, *z, *z, *mean))


# ---------- beginning of script part

# run the computations with multiple cores
ncpus = int(os.environ.get('SLURM_CPUS_PER_TASK', default=1))
print('Beginning computations with {} CPUs...'.format(ncpus))
with multiprocessing.Pool(processes=ncpus) as pool:
    if t is not None:
        stars_results_tuple, formed_results_tuple, gas_results_tuple, halo_results_tuple = zip(
            *pool.starmap(process_snapshot, zip(snaps, t), chunksize=1))
    else:
        stars_results_tuple, formed_results_tuple, gas_results_tuple, halo_results_tuple = zip(
            *pool.map(process_snapshot, snaps, chunksize=1))

# make headers for save files
if t is not None:
    header = delimiter.join((s.ljust(fieldwidth) for s in (
            'Time',
            'Total',
            'N_particles',
            *('s{:.1f}kpc'.format(r) for r in R),
            *('x{:.1f}kpc'.format(r) for r in R),
            *('y{:.1f}kpc'.format(r) for r in R),
            *('z{:.1f}kpc'.format(r) for r in R),
            *('{:.1f}kpc'.format(r) for r in R),
        )))
else:
    header = delimiter.join((s.ljust(fieldwidth) for s in (
            'Total',
            'N_particles',
            *('s{:.1f}kpc'.format(r) for r in R),
            *('x{:.1f}kpc'.format(r) for r in R),
            *('y{:.1f}kpc'.format(r) for r in R),
            *('z{:.1f}kpc'.format(r) for r in R),
            *('{:.1f}kpc'.format(r) for r in R),
        )))

# convert tuples to strings
stars_results_string = header + '\n' + '\n'.join(stars_results_tuple) + '\n'
formed_results_string = header + '\n' + '\n'.join(formed_results_tuple) + '\n'
gas_results_string = header + '\n' + '\n'.join(gas_results_tuple) + '\n'
halo_results_string = header + '\n' + '\n'.join(halo_results_tuple) + '\n'

# write the output
if args.print:
    print('Printing results...\n\n')
    print('===== STELLAR MASS =====')
    print(stars_results_string)
    print('===== STELLAR MASS FORMED =====')
    print(formed_results_string)
    print('===== GAS MASS =====')
    print(gas_results_string)
    print('===== HALO MASS =====')
    print(halo_results_string)
else:
    print('Writing files...')
    with open(args.prefix + '_mass_table_stars.txt' if args.prefix is not None else 'mass_table_stars.txt', 'w') as f:
        f.write(stars_results_string)
    with open(args.prefix + '_mass_table_stars_formed.txt' if args.prefix is not None else 'mass_table_stars_formed.txt', 'w') as f:
        f.write(formed_results_string)
    with open(args.prefix + '_mass_table_gas.txt' if args.prefix is not None else 'mass_table_gas.txt', 'w') as f:
        f.write(gas_results_string)
    with open(args.prefix + '_mass_table_halo.txt' if args.prefix is not None else 'mass_table_halo.txt', 'w') as f:
        f.write(halo_results_string)
    print('Done!')
