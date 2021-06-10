import h5py
import numpy as np
from argparse import ArgumentParser
import os.path

# warning: currently there will be bugs if --mmg is used

parser = ArgumentParser()
parser.add_argument('data_dir')
parser.add_argument('--mmg', action = 'store_true')  # Use most massive galaxy rather than peak brightness
args = parser.parse_args()
data_dir = args.data_dir
if args.mmg:
    particle_ids_center = np.loadtxt(os.path.join(data_dir, 'mmg_particle_ids.txt'))
if args.mmg:
    save_file = os.path.join(data_dir, 'star_mass_data_mmg.txt')
    save_file_x = os.path.join(data_dir, 'star_mass_data_mmg_x.txt')
    save_file_mean = os.path.join(data_dir, 'star_mass_data_mmg_mean.txt')
else:
    save_file = os.path.join(data_dir, 'star_mass_data.txt')
    save_file_x = os.path.join(data_dir, 'star_mass_data_x.txt')
    save_file_mean = os.path.join(data_dir, 'star_mass_data_mean.txt')


# ---------- parameters

snaps = np.arange(0, 3) # numbers of snapshot files to process
R = [30, 50, 70] # radii within which to compute mass
delimiter = '  ' # separator between fields in the output file
precision = 6 # significant figures to use in the output file
fieldwidth = 11 # all fields padded to be at least this width in the output file


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

# compute 2D array of output
output_matrix = np.fromiter(
        generator_fcn(data_dir, snaps),
        dtype='float',
        count=len(snaps)*(2+4*len(R))
    ).reshape((len(snaps), 2+4*len(R)))
# extract only part for cylinder about x-axis
output_matrix_x = output_matrix[:, (0, 1, *range(2+len(R), 2+2*len(R)))]
# compute mean results for the three cylinder orientations
output_matrix_mean = np.concatenate((
        output_matrix[:, (0,1)],
        (output_matrix[:,2+len(R):2+2*len(R)]
         + output_matrix[:,2+2*len(R):2+3*len(R)]
         + output_matrix[:,2+3*len(R):2+4*len(R)])
        / 3.),
    axis=1)
# make headers for save files
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

