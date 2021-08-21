#!/bin/sh

# The purpose of this script is to re-generate the initial conditions in the 
# real*/ directories from the HDF5 files for each galaxy. If we change the
# mass budget and replace the HDF5 file for each galaxy in this directory,
# then we can run this script to re-splice initial conditions while keeping the
# same initial position and velocity for each galaxy.

#module load hdf5 scipy-stack

for real in real01{0..4}; do
    echo "Splicing $real"
    echo "python ../generate_protocluster.py $real.hdf5 --log $real-log --load $real-log --no-copy"
    python ../generate_protocluster.py $real.hdf5 --log $real-log --load $real-log --no-copy
done
