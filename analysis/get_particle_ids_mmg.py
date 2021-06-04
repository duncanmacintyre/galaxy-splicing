import h5py
import numpy as np
import argparse as ap


parser = ap.ArgumentParser()
parser.add_argument('data_dir')
args = parser.parse_args()

# I actually can't grab the ParticleIDs directly from
# A.hdf5 since I reset the ParticleIDs in the splicing
# Python script.
#
# However, I always did the same ordering when I spliced
# together multiple galaxies.
# That ordering is GADGET/GIZMO ordering, so
# PartType0, PartType1, PartType2, PartType3, PartType4
#
# I need the ParticleIDs for PartType2. The first N 
# particles of PartType2 in the SPLICED IC file will
# contain the ParticleIDs for those particles.
#
# In this script, I need to find out how many PartType2
# are in galaxy A so that I can only grab those from the
# initial conditions.
if args.data_dir == 'makegalaxy_g0.5':
    data_dir = 'SPT2349_1e6_gf0.5'
elif args.data_dir == 'makegalaxy_g0.7':
    data_dir = 'SPT2349_1e6_real01'
else:
    data_dir = 'SPT2349_1e6_gf0.9'


with h5py.File('%s/A.hdf5' % args.data_dir, 'r') as f:
    particle_ids = np.array(f['/PartType2/ParticleIDs'], dtype = np.uint32)
    N_stars = len(particle_ids)
    del particle_ids

print('There are %d stars in galaxy A.' % N_stars)

with h5py.File('%s/spliced_1e6.hdf5' % data_dir, 'r') as f:
    particle_ids = np.array(f['/PartType2/ParticleIDs'], dtype = np.uint32)

particle_ids = particle_ids[0:N_stars]

with open('%s_mmg_particle_ids.txt' % args.data_dir, 'w') as f:
    for particle_id in particle_ids:
        f.write('%d\n' % particle_id)

