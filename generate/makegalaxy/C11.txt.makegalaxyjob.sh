#!/bin/bash
#SBATCH --time=00:60:00   # walltime in d-hh:mm or hh:mm:ss format
#SBATCH --job-name="mg C11.txt"
#SBATCH --output=C11.txt.slurm_log
#SBATCH --error=C11.txt.slurm_error
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dm1@student.ubc.ca
#SBATCH --mem=1500  # max memory
echo Time: $(date)
echo "Starting makegalaxy for C11.txt. The job ID is $SLURM_JOBID."
makegalaxy_e1dcb87 C11.txt
