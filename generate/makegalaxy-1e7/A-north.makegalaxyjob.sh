#!/bin/bash
#SBATCH --time=01:10:00   # walltime in d-hh:mm or hh:mm:ss format
#SBATCH --job-name="mg A-north"
#SBATCH --output=A-north.slurm_log
#SBATCH --error=A-north.slurm_error
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dm1@student.ubc.ca
#SBATCH --mem=1200  # max memory
echo Time: $(date)
echo "Starting makegalaxy for A-north. The job ID is $SLURM_JOBID."
echo $SLURM_JOBID >> job_ids
/home/dm1/bin/makegalaxy_ff453fa A-north.txt
