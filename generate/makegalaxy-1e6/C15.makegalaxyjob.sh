#!/bin/bash
#SBATCH --time=01:10:00   # walltime in d-hh:mm or hh:mm:ss format
#SBATCH --job-name="mg C15"
#SBATCH --output=C15.slurm_log
#SBATCH --error=C15.slurm_error
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dm1@student.ubc.ca
#SBATCH --mem=1200  # max memory
echo Time: $(date)
echo "Starting makegalaxy for C15. The job ID is $SLURM_JOBID."
echo $SLURM_JOBID >> job_ids
/home/dm1/bin/makegalaxy_ff453fa C15.txt
