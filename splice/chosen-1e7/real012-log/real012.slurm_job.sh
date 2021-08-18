#!/bin/bash
#SBATCH --time=1-00:00   # walltime in d-hh:mm or hh:mm:ss format
#SBATCH --nodes=4
#SBATCH --mem=0
#SBATCH --ntasks-per-node=48
#SBATCH --job-name="real012"
#SBATCH --output="temp_slurm_log"
#SBATCH --error="temp_slurm_error"
#SBATCH --mail-type=START,END,FAIL
#SBATCH --mail-user={email}
#####       SBATCH --constraint=[dragonfly1|dragonfly2|dragonfly3|dragonfly4|dragonfly5]     # only turn on this line on Niagara

# this variable is 1 for first time Gizmo run, 2 for first resubmission, 3 for second resubmission, ...
declare -i this_job_number

# the how-many-jobs-run file will contain the number of GIZMO jobs that have ran for this simulation
# has this simulation previously been running? to decide, look for the how-many-jobs-run file
if test -f how-many-jobs-run; then # we have run before
    this_job_number=$(cat how-many-jobs-run)
    this_job_number+=1
else # this is the first run
    this_job_number=1

    # let's set up the resubmit script that Gizmo will call
    echo "#!/bin/sh" > resubmit.sh
    echo "touch resubmit-flag" >> resubmit.sh
    chmod u+x resubmit.sh
fi

echo $this_job_number > how-many-jobs-run
echo $SLURM_JOBID >> job-ids

# add header markers to the log files
echo "" >> gizmo-log
echo "======================= real012 Slurm run $this_job_number, job ID $SLURM_JOBID =======================" >> gizmo-log
echo "" >> gizmo-error
echo "======================= real012 Slurm run $this_job_number, job ID $SLURM_JOBID =======================" >> gizmo-error
echo "" >> slurm-log
echo "======================= real012 Slurm run $this_job_number, job ID $SLURM_JOBID =======================" >> slurm-log
echo "" >> slurm-error
echo "======================= real012 Slurm run $this_job_number, job ID $SLURM_JOBID =======================" >> slurm-error

# start GIZMO
echo Time: $(date)
echo "Starting GIZMO for real012 on $SLURM_JOB_NUM_NODES nodes ($SLURM_NTASKS tasks). The job ID is $SLURM_JOBID."
echo "Node list: $SLURM_JOB_NODELIST"
#module load CCEnv arch/avx2 nixpkgs/16.09 intel/2016.4 openmpi/2.1.1 fftw-mpi/2.1.5 grackle/3.1 gsl/2.2.1 # on Niagara
module load nixpkgs/16.09 intel/2016.4 openmpi/2.1.1 fftw-mpi/2.1.5 grackle/3.1 gsl/2.2.1
echo "Loaded modules:"
module list
echo ""
if test $this_job_number -eq 1; then # starting for first time
    echo "Starting GIZMO for the first time."
    mpirun $HOME/bin/GIZMO_d768f3d real012.param >> gizmo-log 2>> gizmo-error
else # start from restart files
    echo "Starting GIZMO from restart files."
    mpirun $HOME/bin/GIZMO_d768f3d real012.param 1 >> gizmo-log 2>> gizmo-error
fi
echo "GIZMO done!"

# submit a new batch job if the resubmit-flag file is present
if test -f resubmit-flag; then
    echo "Resubmitting the simulation to Slurm."
    sbatch real012.slurm_job.sh
    echo "Removing the resubmission flag."
    rm resubmit-flag
fi

# append this run's slurm log and error to main one
cat temp_slurm_log >> slurm-log
cat temp_slurm_error >> slurm-error

echo "Done!"
