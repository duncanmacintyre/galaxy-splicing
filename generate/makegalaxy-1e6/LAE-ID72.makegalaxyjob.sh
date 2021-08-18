#!/bin/bash
#SBATCH --time=01:10:00   # walltime in d-hh:mm or hh:mm:ss format
#SBATCH --job-name="mg LAE-ID72"
#SBATCH --output=LAE-ID72.slurm_log
#SBATCH --error=LAE-ID72.slurm_error
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dm1@student.ubc.ca
#SBATCH --mem=1200  # max memory
echo Time: $(date)
echo "Starting makegalaxy for LAE-ID72. The job ID is $SLURM_JOBID."
echo $SLURM_JOBID >> job_ids
/home/dm1/bin/makegalaxy_ff453fa LAE-ID72.txt
