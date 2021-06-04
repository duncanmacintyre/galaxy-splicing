import h5py
import numpy as np
import argparse as ap


parser = ap.ArgumentParser()
parser.add_argument('data_dir')
parser.add_argument('--mmg', action = 'store_true')  # Use most massive galaxy rather than peak brightness
args = parser.parse_args()

if args.data_dir == 'makegalaxy_g0.5':
    data_dir = 'SPT2349_1e6_gf0.5'
elif args.data_dir == 'makegalaxy_g0.7':
    data_dir = 'SPT2349_1e6_real01'
else:
    data_dir = 'SPT2349_1e6_gf0.9'

snaps = np.arange(0, 103)

total_formed_mass = np.zeros(len(snaps))
total_stellar_mass = np.copy(total_formed_mass)
mass_in_5kpc = np.copy(total_formed_mass)
mass_in_7kpc = np.copy(total_formed_mass)
mass_in_15kpc = np.copy(total_formed_mass)
mass_in_30kpc = np.copy(total_formed_mass)
mass_in_50kpc = np.copy(total_formed_mass)
mass_in_70kpc = np.copy(total_formed_mass)

mass_in_5kpc_b = np.copy(total_formed_mass)
mass_in_7kpc_b = np.copy(total_formed_mass)
mass_in_15kpc_b = np.copy(total_formed_mass)
mass_in_30kpc_b = np.copy(total_formed_mass)
mass_in_50kpc_b = np.copy(total_formed_mass)
mass_in_70kpc_b = np.copy(total_formed_mass)

mass_in_5kpc_c = np.copy(total_formed_mass)
mass_in_7kpc_c = np.copy(total_formed_mass)
mass_in_15kpc_c = np.copy(total_formed_mass)
mass_in_30kpc_c = np.copy(total_formed_mass)
mass_in_50kpc_c = np.copy(total_formed_mass)
mass_in_70kpc_c = np.copy(total_formed_mass)

if args.mmg:
    particle_ids_center = np.loadtxt('%s_mmg_particle_ids.txt' % args.data_dir)

def get_mass_in_cylinder(masses, coords, column_to_remove):
    radii = np.linalg.norm(np.delete(coords, column_to_remove, 1), axis = 1)

    # 5 kpc radius
    in_region = np.where(radii < 5)
    mass_in_5 = np.sum(masses[in_region])

    # 7 kpc radius
    in_region = np.where(radii < 7)
    mass_in_7 = np.sum(masses[in_region])

    # 15 kpc radius
    in_region = np.where(radii < 15)
    mass_in_15 = np.sum(masses[in_region])

    # 30 kpc radius
    in_region = np.where(radii < 30)
    mass_in_30 = np.sum(masses[in_region])

    # 50 kpc radius
    in_region = np.where(radii < 50)
    mass_in_50 = np.sum(masses[in_region])

    # 70 kpc radius
    in_region = np.where(radii < 70)
    mass_in_70 = np.sum(masses[in_region])

    if column_to_remove == 2:
        mass_in_5kpc[i] = mass_in_5
        mass_in_7kpc[i] = mass_in_7
        mass_in_15kpc[i] = mass_in_15
        mass_in_30kpc[i] = mass_in_30
        mass_in_50kpc[i] = mass_in_50
        mass_in_70kpc[i] = mass_in_70
    elif column_to_remove == 1:
        mass_in_5kpc_b[i] = mass_in_5
        mass_in_7kpc_b[i] = mass_in_7
        mass_in_15kpc_b[i] = mass_in_15
        mass_in_30kpc_b[i] = mass_in_30
        mass_in_50kpc_b[i] = mass_in_50
        mass_in_70kpc_b[i] = mass_in_70
    else:
        mass_in_5kpc_c[i] = mass_in_5
        mass_in_7kpc_c[i] = mass_in_7
        mass_in_15kpc_c[i] = mass_in_15
        mass_in_30kpc_c[i] = mass_in_30
        mass_in_50kpc_c[i] = mass_in_50
        mass_in_70kpc_c[i] = mass_in_70

for i, snap in enumerate(snaps):
    snap_str = str(snap).zfill(3)

    data_file = '%s/snapshot_%s.hdf5' % (data_dir, snap_str)

    print('Operating on %s' % data_file)

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
            mmg_idx = np.in1d(disk_ids, particle_ids_center)
            new_mmg_idx = np.array([not i for i in mmg_idx])
            mmg_idx = new_mmg_idx
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
    
    get_mass_in_cylinder(masses, coords, 0)  # yz plane
    get_mass_in_cylinder(masses, coords, 1)  # xz plane
    get_mass_in_cylinder(masses, coords, 2)  # xy plane

    total_formed_mass[i] = np.sum(stellar_masses)
    total_stellar_mass[i] = np.sum(disk_masses) + total_formed_mass[i]

if args.mmg:
    save_file = '%s_star_mass_data_mmg.txt' % data_dir
else:
    save_file = '%s_star_mass_data.txt' % data_dir

with open(save_file, 'w') as f:
    f.write('Formed\tTotal\ta5kpc\ta7kpc\ta15kpc\a30kpc\ta50kpc\ta70kpc\tb5kpc\tb7kpc\tb15kpc\tb30kpc\tb50kpc\tb70kpc\tc5kpc\tc7kpc\tc15kpc\tc30kpc\tc50kpc\tc70kpc\n')
    for i in range(len(snaps)):
        f.write('%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\n' % (total_formed_mass[i], 
                                                                  total_stellar_mass[i],
                                                                  mass_in_5kpc[i], mass_in_7kpc[i], mass_in_15kpc[i],
                                                                  mass_in_30kpc[i], mass_in_50kpc[i], mass_in_70kpc[i],
                                                                  mass_in_5kpc_b[i], mass_in_7kpc_b[i], mass_in_15kpc_b[i],
                                                                  mass_in_30kpc_b[i], mass_in_50kpc_b[i], mass_in_70kpc_b[i],
                                                                  mass_in_5kpc_c[i], mass_in_7kpc_c[i], mass_in_15kpc_c[i],
                                                                  mass_in_30kpc_c[i], mass_in_50kpc_c[i], mass_in_70kpc_c[i]))
