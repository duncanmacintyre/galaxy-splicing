#!/bin/bash
#SBATCH --time=00:40:00   # walltime in d-hh:mm or hh:mm:ss format
#SBATCH --job-name="mg test/C7.txt"
#SBATCH --mem=2G    # max memory
./makegalaxy test/C7.txt
