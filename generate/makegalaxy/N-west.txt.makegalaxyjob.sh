#!/bin/bash
#SBATCH --time=00:40:00   # walltime in d-hh:mm or hh:mm:ss format
#SBATCH --job-name="mg N-west.txt"
#SBATCH --output=N-west.txt.slurm_log
#SBATCH --error=N-west.txt.slurm_error
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dm1@student.ubc.ca
#SBATCH --mem=1300  # max memory
echo Time: $(date)
echo "Starting makegalaxy for N-west.txt. The job ID is $SLURM_JOBID."
makegalaxy_e1dcb87 N-west.txt
