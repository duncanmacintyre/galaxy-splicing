#!/bin/bash
#SBATCH --time=00:40:00   # walltime in d-hh:mm or hh:mm:ss format
#SBATCH --job-name="mg C7.txt"
#SBATCH --output=C7.txt.slurm_log
#SBATCH --error=C7.txt.slurm_error
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dm1@student.ubc.ca
#SBATCH --mem=1300  # max memory
echo Time: $(date)
echo "Starting makegalaxy for C7.txt. The job ID is $SLURM_JOBID."
makegalaxy_e1dcb87 C7.txt
