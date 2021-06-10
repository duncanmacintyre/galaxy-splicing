#!/bin/python
# DR: original version of code (for 2020 Rennehan et al. paper)
# DM: rewrote to make it easy to change radii, added more output files, refactored
#     for speed, readability, and versatility (June 2021)
#
# use --help for usage
# warning: currently there will be bugs if --mmg is used
#

import h5py
import numpy as np
from argparse import ArgumentParser
import os.path
import sys

# ---------- argument parsing

parser = ArgumentParser()
parser.add_argument('data_dir', help='path to directory containing snapshot files')
g1 = parser.add_mutually_exclusive_group(required=True)
g1.add_argument('-R', nargs='+', type=float, help='radii to use')
g1.add_argument('-r', nargs=2, metavar=('DR', 'MAX_R'), type=float, help='will use radii DR, 2DR, 3DR, ... up to but not including MAX_R')
g2 = parser.add_mutually_exclusive_group(required=True)
g2.add_argument('-N', nargs='+', type=int, help='the snapshot file numbers to use')
g2.add_argument('-n', type=int, metavar='MAX_SNAP', help='use snapshot numbers 0, 1, 2, 3, ... up to but not including MAX_SNAP; e.g. use 100 for 0 through 99')
g3 = parser.add_mutually_exclusive_group(required=False)
g3.add_argument('-T', nargs='+', type=float, help='the times for each snapshot number')
g3.add_argument('-t', metavar='DT', type=float, help='snapshots have time N*DT where N is snapshot number')
parser.add_argument('--mmg', action = 'store_true', help='use most massive galaxy for centre rather than peak brightness; not currently working')
parser.add_argument('-d', default='  ', help='separator between fields in the output files; default two spaces')
parser.add_argument('-p', default=6, type=int, help='how many significant figures to use in the output file; default 6')
parser.add_argument('-w', default=11, type=int, help='all fields padded to be at least this width in the output file; default 11')
args = parser.parse_args()

# data_dir - the path to the folder containing snapshot files
data_dir = args.data_dir

# R - the radii we will use
if args.R is not None: # -R was given
    R = np.array(args.R)
else: # -r was given
    R = np.arange(args.r[0], args.r[1], args.r[0])
assert(len(R) > 0)
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
    t = np.array(args.T).reshape(-1, 1)
    assert(all(t >= 0))
    assert(len(t) == len(snaps))
elif args.t is not None:
    assert(args.t > 0)
    t = (args.t * snaps).reshape(-1, 1)
else:
    t = None

# initialize variables based on mmg
if args.mmg:
    particle_ids_center = np.loadtxt(os.path.join(data_dir, 'mmg_particle_ids.txt'))
    save_file = os.path.join(data_dir, 'star_mass_data_mmg.txt')
    save_file_x = os.path.join(data_dir, 'star_mass_data_mmg_x.txt')
    save_file_mean = os.path.join(data_dir, 'star_mass_data_mmg_mean.txt')
else:
    save_file = os.path.join(data_dir, 'star_mass_data.txt')
    save_file_x = os.path.join(data_dir, 'star_mass_data_x.txt')
    save_file_mean = os.path.join(data_dir, 'star_mass_data_mean.txt')

# output formatting
delimiter = args.d # separator between fields in the output file
precision = args.p # significant figures to use in the output file
assert(precision > 0)
fieldwidth = args.w # all fields padded to be at least this width in the output file
assert(fieldwidth > 0)


# ---------- functions

# given particle masses and coords, compute stellar mass contained within radius r
#   column_to_remove=None         sphere
#   column_to_remove=0            cylinder about x-axis
#   column_to_remove=1            cylinder about y-axis
#   column_to_remove=2            cylinder about z-axis
def get_mass_in_cylinder(masses, coords, r, column_to_remove):
    if column_to_remove is None:
        radii = np.linalg.norm(coords, axis=1)
    else:
        radii = np.linalg.norm(np.delete(coords, column_to_remove, axis=1), axis=1)

    return(np.sum(masses[radii < r]))


# GENERATOR: yield output for all snapshots
# yields 2+4*len(R) consecutive items for each snapshot, with...
#   first item: stellar mass formed so far
#   second item: total stellar mass
#   next len(R) items: mass within spheres of radii from R
#   next len(R) items: mass within cylinders about x-axis of radii from R
#   next len(R) items: mass within cylinders about y-axis of radii from R
#   next len(R) items: mass within cylinders about z-axis of radii from R
# in total, yields N*(2+4*len(R)) items, where N is the number of snapshots
#
# dirname is the path to the folder containing the snapshots
# snaps is an iterator giving snapshot numbers, e.g. use range(100) for snapshots 0 through 99
def generator_fcn(data_dir, snaps):
    for snap in snaps:

        data_file = os.path.join(data_dir, 'snapshot_{}.hdf5'.format(str(snap).zfill(3)))
        print('Operating on {}'.format(data_file))
        sys.stdout.flush()

        with h5py.File(data_file, 'r') as f:
            disk_ids = np.array(f['/PartType2/ParticleIDs'])
            disk_masses = np.array(f['/PartType2/Masses']) * 1e10
            disk_coords = np.array(f['/PartType2/Coordinates'])
            try:
                stellar_masses = np.array(f['/PartType4/Masses']) * 1e10
                stellar_coords = np.array(f['/PartType4/Coordinates'])
            except:
                stellar_masses = np.zeros(len(disk_masses))
                stellar_coords = None

        if stellar_coords is not None: 
            coords = np.vstack((disk_coords, stellar_coords))
        else:
            coords = disk_coords

        if args.mmg:
            if data_dir == 'SPT2349_1e6_gf0.9':
                mmg_idx = np.array([not i for i in np.in1d(disk_ids, particle_ids_center)])
            else:
                mmg_idx = np.in1d(disk_ids, particle_ids_center)
            mmg_disk_masses = disk_masses[mmg_idx]
            mmg_disk_coords = disk_coords[mmg_idx]

            mmg_mass_sum = np.sum(mmg_disk_masses)
            mmg_com_x = np.sum(mmg_disk_masses * mmg_disk_coords[:, 0]) / mmg_mass_sum
            mmg_com_y = np.sum(mmg_disk_masses * mmg_disk_coords[:, 1]) / mmg_mass_sum
            mmg_com_z = np.sum(mmg_disk_masses * mmg_disk_coords[:, 2]) / mmg_mass_sum
            
            center_offset = np.array([mmg_com_x, mmg_com_y, mmg_com_z])
        else:
            bins = 512
            low_limit = -75
            up_limit = 75

            bounds = [[low_limit, up_limit], [low_limit, up_limit], [low_limit, up_limit]]
            map_slope = (up_limit - low_limit) / bins
            map_offset = low_limit

            hist, be = np.histogramdd(coords, bins = bins, range = bounds)
            pos = np.argwhere(hist == hist.max())[0]
            ruler = np.linspace(low_limit, up_limit, bins)
            pos_x = map_slope * pos[0] + map_offset
            pos_y = map_slope * pos[1] + map_offset
            pos_z = map_slope * pos[2] + map_offset

            center_offset = np.array([pos_x, pos_y, pos_z])

        # Center on the peak brightness or most massive galaxy
        coords -= center_offset

        if stellar_coords is not None:
            masses = np.concatenate((disk_masses, stellar_masses))
        else:
            masses = disk_masses

        yield(np.sum(stellar_masses)) # formed stellar mass
        yield(np.sum(disk_masses) + np.sum(stellar_masses)) # total stellar mass
        yield from (get_mass_in_cylinder(masses, coords, r, None) for r in R) # masses within spheres
        yield from (get_mass_in_cylinder(masses, coords, r, 0) for r in R) # masses within spheres cylinder about x-axis
        yield from (get_mass_in_cylinder(masses, coords, r, 1) for r in R) # masses within spheres cylinder about y-axis
        yield from (get_mass_in_cylinder(masses, coords, r, 2) for r in R) # masses within spheres cylinder about z-axis


# ---------- beginning of script part

print('Beginning computations...')

# number of columns in output before the radial bin columns
n_head_cols = 3 if t is not None else 2

# compute 2D array of output
output_matrix = np.fromiter(
        generator_fcn(data_dir, snaps),
        dtype='float',
        count=len(snaps)*(2+4*len(R))
    ).reshape((len(snaps), 2+4*len(R)))

# if time is set, add a column for time
if t is not None:
    output_matrix = np.concatenate((t, output_matrix), axis=1)

# extract only part for cylinder about x-axis
output_matrix_x = output_matrix[:, (*range(n_head_cols), *range(n_head_cols+len(R), n_head_cols+2*len(R)))]
# compute mean results for the three cylinder orientations
output_matrix_mean = np.concatenate((
        output_matrix[:, 0:n_head_cols],
        (output_matrix[:,n_head_cols+len(R):n_head_cols+2*len(R)]
         + output_matrix[:,n_head_cols+2*len(R):n_head_cols+3*len(R)]
         + output_matrix[:,n_head_cols+3*len(R):n_head_cols+4*len(R)])
        / 3.),
    axis=1)

# make headers for save files
if t is not None:
    header = delimiter.join((s.ljust(fieldwidth) for s in (
            'Time',
            'Formed',
            'Total',
            *("s{:.1f}kpc".format(r) for r in R),
            *("x{:.1f}kpc".format(r) for r in R),
            *("y{:.1f}kpc".format(r) for r in R),
            *("z{:.1f}kpc".format(r) for r in R),
        )))
    header_x = delimiter.join((s.ljust(fieldwidth) for s in (
            'Time',
            'Formed',
            'Total',
            *("{:.1f}kpc".format(r) for r in R),
        )))
else:
    header = delimiter.join((s.ljust(fieldwidth) for s in (
            'Formed',
            'Total',
            *("s{:.1f}kpc".format(r) for r in R),
            *("x{:.1f}kpc".format(r) for r in R),
            *("y{:.1f}kpc".format(r) for r in R),
            *("z{:.1f}kpc".format(r) for r in R),
        )))
    header_x = delimiter.join((s.ljust(fieldwidth) for s in (
            'Formed',
            'Total',
            *("{:.1f}kpc".format(r) for r in R),
        )))
header_mean = header_x

# save to txt files
np.savetxt(save_file, output_matrix,
           fmt='%{:d}.{:d}g'.format(fieldwidth, precision), delimiter=delimiter,
           header=header, comments='')
np.savetxt(save_file_x, output_matrix_x,
           fmt='%{:d}.{:d}g'.format(fieldwidth, precision), delimiter=delimiter,
           header=header_x, comments='')
np.savetxt(save_file_mean, output_matrix_mean,
           fmt='%{:d}.{:d}g'.format(fieldwidth, precision), delimiter=delimiter,
           header=header_mean, comments='')

