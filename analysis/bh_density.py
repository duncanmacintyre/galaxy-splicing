import yt
import numpy as np
import matplotlib.pyplot as plt

def snap2time(snapnum):
    unit_time = 978.028 # In Myr
    dt = 0.01 # Time between snapshots
    myr_to_gyr = 1.0e-3
    
    return snapnum * dt * unit_time * myr_to_gyr

num_snaps = 103

bbox_lim = 1e4
bbox_vec = [-bbox_lim, bbox_lim]
bbox = [bbox_vec, bbox_vec, bbox_vec]

n_bhs_2kpc = np.zeros(num_snaps, dtype = int)
n_bhs_1kpc = np.copy(n_bhs_2kpc)

mass_bhs_2kpc = np.zeros(num_snaps, dtype = float)
mass_bhs_1kpc = np.copy(mass_bhs_2kpc)

for i in range(1, num_snaps):
    snap = str(i).zfill(3)
    data_file = 'SPT2349_1e6_real01/snapshot_%s.hdf5' % snap
    
    ds = yt.load(data_file, bounding_box = bbox)
    ad = ds.all_data()
   
    star_coords = ad[('PartType4', 'Coordinates')].in_units('kpc')
    star_masses = ad[('PartType4', 'Masses')].in_units('Msun')
    bh_coords = ad[('PartType5', 'Coordinates')].in_units('kpc')
    bh_masses = ad[('PartType5', 'Masses')].in_units('Msun')
    
    com_x = np.sum(star_masses * star_coords[:, 0]) / np.sum(star_masses)
    com_y = np.sum(star_masses * star_coords[:, 1]) / np.sum(star_masses)
    com_z = np.sum(star_masses * star_coords[:, 2]) / np.sum(star_masses)
    
    star_coords[:, 0] -= com_x
    star_coords[:, 1] -= com_y
    star_coords[:, 2] -= com_z
    
    bh_coords[:, 0] -= com_x
    bh_coords[:, 1] -= com_y
    bh_coords[:, 2] -= com_z
    
    star_distances = np.linalg.norm(star_coords, axis = 1)
    bh_distances = np.linalg.norm(bh_coords, axis = 1)
    
    bh_in_2kpc_idx = np.where(bh_distances < 2)
    bh_in_1kpc_idx = np.where(bh_distances < 1)
    
    n_bhs_2kpc[i] = len(bh_in_2kpc_idx[0])
    n_bhs_1kpc[i] = len(bh_in_1kpc_idx[0])
    
    mass_bhs_2kpc[i] = np.sum(bh_masses[bh_in_2kpc_idx])
    mass_bhs_1kpc[i] = np.sum(bh_masses[bh_in_1kpc_idx])

plt.figure()
plt.xlim([0, 1])
plt.ylim([0, 14])
plt.xlabel('Time [Gyr]')
plt.ylabel('Total Blackhole Count')
plt.plot(snap2time(np.arange(0, num_snaps)), n_bhs_2kpc, label = '<2kpc')
plt.plot(snap2time(np.arange(0, num_snaps)), n_bhs_1kpc, label = '<1kpc')
plt.legend(loc = 'upper left')
plt.savefig('bh_counts.png')
plt.close()

plt.figure()
plt.xlim([0, 1])
plt.ylim([1e5, 1e8])
plt.yscale('log')
plt.xlabel('Time [Gyr]')
plt.ylabel('Total Blackhole Mass [Msun]')
plt.plot(snap2time(np.arange(0, num_snaps)), mass_bhs_2kpc, label = '<2kpc')
plt.plot(snap2time(np.arange(0, num_snaps)), mass_bhs_1kpc, label = '<1kpc')
plt.legend(loc = 'upper left')
plt.savefig('bh_masses.png')
plt.close()
