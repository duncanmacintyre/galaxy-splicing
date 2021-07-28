#!/bin/sh
#SBATCH --ntasks=1
#SBATCH --time=04:00:00
#SBATCH --mem=4000
#SBATCH --out=make_some_ics_log
splice 001.hdf5 $(cat galaxies) --no-copy
splice 002.hdf5 $(cat galaxies) --no-copy
splice 003.hdf5 $(cat galaxies) --no-copy
splice 004.hdf5 $(cat galaxies) --no-copy
splice 005.hdf5 $(cat galaxies) --no-copy
splice 006.hdf5 $(cat galaxies) --no-copy
splice 007.hdf5 $(cat galaxies) --no-copy
splice 008.hdf5 $(cat galaxies) --no-copy
splice 009.hdf5 $(cat galaxies) --no-copy
splice 010.hdf5 $(cat galaxies) --no-copy
splice 011.hdf5 $(cat galaxies) --no-copy
splice 012.hdf5 $(cat galaxies) --no-copy
splice 013.hdf5 $(cat galaxies) --no-copy
splice 014.hdf5 $(cat galaxies) --no-copy
splice 015.hdf5 $(cat galaxies) --no-copy
