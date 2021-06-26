
import os
import sys
import math
import h5py

gizmo_top_dir = os.path.abspath('/scratch/{}/gizmo-galaxies-in-isolation'.format(os.getenv('USER')))
n_tasks = 4

# ----- templates -----

template_slurm_script = """#!/bin/bash
#SBATCH --time={3}   # walltime in d-hh:mm or hh:mm:ss format
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={4}
#SBATCH --mem={2}  # max memory
#SBATCH --job-name="gz {0}"
#SBATCH --output={0}.slurm_log
#SBATCH --error={0}.slurm_error
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=dm1@student.ubc.ca
echo Time: $(date)
echo "Starting makegalaxy for {0}. The job ID is $SLURM_JOBID."
echo $SLURM_JOBID >> job_ids
{1} {0}.param

"""

template_param_file = """
%-------------------------------------------------------------------------
%----  This file contains the input parameters needed at run-time for 
%       simulations. It is based on and closely resembles the GADGET-3
%       parameterfile (format of which and parsing routines written by 
%       Volker Springel [volker.springel@h-its.org]). It has been updated
%       with new naming conventions and additional variables as needed by 
%       Phil Hopkins [phopkins@caltech.edu] for GIZMO.
%-------------------------------------------------------------------------

%----  Relevant files
InitCondFile {0}
%InitCondFile ./snapshot_001
OutputDir     ./

%---- File formats 
ICFormat    3  % 1=binary, 3=hdf5, 4=cluster 
SnapFormat  3  % 1=binary, 3=hdf5 

%---- Output parameters 
RestartFile                 restart
SnapshotFileBase            snapshot
OutputListOn                0  % =1 to use list in "OutputListFilename" 
OutputListFilename          output_times.txt  % list of times for snaps 
NumFilesPerSnapshot         1
NumFilesWrittenInParallel   1  % must be < N_processors & power of 2

%---- Output frequency 
TimeOfFirstSnapshot     0.0
TimeBetSnapshot         0.01 
TimeBetStatistics       0.01 

%---- CPU-time limits 
TimeLimitCPU            {2}   %  259200  % in seconds 
CpuTimeBetRestartFile   200    %  3500    % in seconds 
ResubmitOn              0
ResubmitCommand         job.bash

%----- Memory alloction 
MaxMemSize          {1}     % sets maximum MPI process memory use in MByte 
PartAllocFactor     3.0     % memory load allowed for better cpu balance 
BufferSize          100      % in MByte 

%---- Characteristics of run 
TimeBegin   0.0    % Beginning of the simulation 
TimeMax     0.25   % End of the simulation 

%---- Cosmological parameters 
ComovingIntegrationOn   0       % is it cosmological? (yes=1, no=0)
BoxSize                 10000.  % in code units
Omega0                  0       % =0 for non-cosmological
OmegaLambda             0       % =0 for non-cosmological
OmegaBaryon             0       % =0 for non-cosmological
HubbleParam             1.0     % little 'h'; =1 for non-cosmological runs

%---- Accuracy of time integration 
MaxSizeTimestep         5.0e-3   % in code units, set for your problem
MinSizeTimestep         1.0e-12 % set this very low, or risk stability

%---- Tree algorithm, force accuracy, domain update frequency 
TreeDomainUpdateFrequency   0.05        % 0.0005-0.05, dept on core+particle number  

%---- System of units 
UnitLength_in_cm            3.085678e21     % 1.0 kpc/h
UnitMass_in_g               1.989e43        % 1.0e10 solar masses/h
UnitVelocity_in_cm_per_s    1.0e5           % 1 km/sec
UnitMagneticField_in_gauss  1.0             % 1 gauss
GravityConstantInternal     0               % calculated by code if =0

%---- Initial temperature & temperature floor 
InitGasTemp     0.      % set by IC file if =0 
MinGasTemp      10.     % don't set <10 in explicit feedback runs, otherwise 0

%---- Density/volume estimation (kernel) 
DesNumNgb               32      % 32 for standard kernel, 60-114 for quintic 
MaxHsml                 1000  % some very large value (here to prevent errors)
MinHsml                 0       % minimum gas kernel length  (=0, should be <=SofteningGasMaxPhys)
MinGasHsmlFractional    1.0

SofteningGas    0.5    % gas (type=0) (in units above, =1 pc softening)
SofteningHalo   0.5 % dark matter/collisionless particles (type=1)
SofteningDisk   1.0   % collisionless particles (type=2)
SofteningBulge  1.0  % collisionless particles (type=3)
SofteningStars  0.5  % stars spawned from gas (type=4)
SofteningBndry  0.5  % black holes (if active), or collisionless (type=5)
%---- if these are set in cosmo runs, SofteningX switches from comoving to physical
%------- units when the comoving value exceeds the choice here
SofteningGasMaxPhys     0.5    % switch to 0.5pc physical below z=1 
SofteningHaloMaxPhys    0.5
SofteningDiskMaxPhys    1.0
SofteningBulgeMaxPhys   1.0
SofteningStarsMaxPhys   0.5
SofteningBndryMaxPhys   0.5

%----- parameters for adaptive gravitational softening 
AGS_DesNumNgb           32  % neighbor number for calculating adaptive gravsoft

%---- Isolated system parameters (for if ISOLATED_SYSTEM is on)
IsolatedSimBlackHoleStart    0.1       % start time for BH feedback, probably ~0.5 is fine
IsolatedSimBlackHoleMass     0.005     % mass of one black hole particle (code units)
IsolatedSimGasFraction       {3}       % initial gas mass fraction (gas mass / (gas mass + star mass))
IsolatedSimMstar             {4}       % total initial stellar mass (in Msun)

%---- Star Formation parameters (GALSF on)
CritPhysDensity     0.2    %  critical physical density for star formation (cm^(-3)) 
SfEffPerFreeFall    0.02   %  SFR/(Mgas/tfreefall) for gas which meets SF criteria 
InitMetallicity     0.01 % Fraction of solar metallicity

%---- sub-grid (Springel+Hernquist/GADGET/AREPO) "effective equation of state" 
%------- star formation+feedback model (GALSF_EFFECTIVE_EQS on)
MaxSfrTimescale     2.8         % code units (SF timescale at 2-phase threshold)
TempSupernova       4.87e7       % in Kelvin (temp of hot gas in 2-phase model)
TempClouds          10000.0      % in Kelvin (temp of cold gas in 2-phase model)
FactorSN            0.18         % SNe coupling frac (frac of egy retained in hot)
FactorEVP           3000.0      % controls Kennicutt normalization 
FactorForSofterEQS  1.0         % interpolate between 'stiff' and isothermal EOS
%------- the sub-grid "decoupled winds" model (GALSF_SUBGRID_WINDS on)
WindEfficiency          2.0     % Not used with Illustris model 
WindEnergyFraction      0.0519  % Value used in old Illustris, a bit smaller than TNG
WindFreeTravelMaxTime   0.1     % 'free-stream time' in units of t_Hubble(z)
WindFreeTravelDensFac   0.01     % 'free-stream' until density < this * CritPhysDensity

%-------------- Black Hole accretion & formation (BLACK_HOLES on)
TimeBetOnTheFlyFoF           1.1            % time (in sec) between FoF searches --> DAA: this is t  % time (in sec) between FoF searches --> DAA: this is t  %%% probably needs to be made small! not in seconds?
BlackHoleAccretionFactor     1.0            % multiplier for mdot
BlackHoleEddingtonFactor     1.0           % fraction of eddington to cap (can be >1)
BlackHoleNgbFactor           4              % multiplier for kernel neighbors for BH
BlackHoleMaxAccretionRadius  2              % max radius for BH neighbor search/accretion
BlackHoleRadiativeEfficiency 0.1        % radiative efficiency
BlackHoleFeedbackFactor      0.05           % generic feedback strength multiplier 
SeedBlackHoleMass            0
BAL_v_outflow                1e4        % v_wind in km/s
VariableWindVelFactor        5.0        % Zhu & Li 2016 Apj (MFM comparison)
VariableWindSpecMomentum     0.

%-------------- Grackle UVB file (GRACKLE on)
GrackleDataFile                 CloudyData_UVB=HM2012_shielded.h5



%-------------- DM - missing fields added
%BAL_f_accretion              0.1        % fraction of gas swallowed by BH (BH_BAL_WINDS) 
%SeedBlackHoleMassSigma       0
%SeedBlackHoleMinRedshift     0
MinFoFMassForNewSeed          100 
%massDMpart                   0.0126849

"""

# ----- functions -----

# given a .hdf5 output file from makegalaxy, print particle counts and masses
#
# arguments:
#       * mass of one disk particle (in Msun)
#       * path to .hdf5 file to analyze
#       * number of tasks
#
# returns tuple of six strings:
#       * gas fraction
#       * total disk mass (in Msun)
#       * suggested memory per task based on number of particles (in MB)
#       * total memory (above times number tasks, in MB)
#       * suggested time limit in seconds
#       * suggested time limit in hh:mm:ss form
#
# We assume (1) no bulge particles, (2) gas and disk particles have same mass per particle.
#
def get_mem_time_info(mass_per_disk_particle_in_Msun, fname, p):

    with h5py.File(fname, 'r') as f:
        N_gas = f['/PartType0/Coordinates'].shape[0]
        N_disk = f['/PartType2/Coordinates'].shape[0]
        N_halo = f['/PartType1/Coordinates'].shape[0]

    # ---- gas fraction ----
    gas_fraction = '{:.8f}'.format(N_gas / (N_gas + N_disk))

    # ---- total disk mass (in Msun) ----
    disk_mass = '{:.9g}'.format(N_disk * mass_per_disk_particle_in_Msun)

    # ---- suggested memory per task based on number of particles (in MB) ---
    m = 300 # horizontal scaling factor
    M = 2000 # max possible memory value we might assign per task (in MB)
    x = N_gas + N_disk + N_halo # total number particles
    memory = (
                 ((M-500)/1500)
                 * (
                     math.log10((m*x/(p*(M-500)) + 400)/100000)
                     / math.log10((m*x/(p*(M-500)) + 400)/100000 + 1)
                     + 1499
                   )
                 + 500
             )
    memory_string = '{:.0f}'.format(memory+2) # per task
    memory_total = '{:.0f}'.format(memory*p) # total 

    # ---- suggested time limit ----
    time_limit_in_s = int(60*15 + 0.1 * x / math.sqrt(p)) # in s
    time_string = str(time_limit_in_s)
    # from https://stackoverflow.com/a/775075/13326516
    m, s = divmod(time_limit_in_s, 60)
    h, m = divmod(m, 60)
    time_hms = '{:02d}:{:02d}:{:02d}'.format(h, m, s) # in hh:mm:ss form

    return (gas_fraction, disk_mass, memory_string, memory_total, time_string, time_hms)




# ----- beginning of script part -----

# find file to which the symlink /home/dm1/bin/GIZMO points
# (this destination will include the commit number)
exec_path = os.readlink('/home/dm1/bin/GIZMO')

# iterate over .hdf5 files; fname is the file name without path or extension
for fname in (os.path.splitext(os.path.basename(s))[0] for s in sys.argv[1:]):

    # path to directory in which we'll store outputs
    dir_path = os.path.join(gizmo_top_dir, fname)

    if os.path.exists(dir_path):
        print('ERROR: {} already exists, skipping {}'.format(dir_path, fname))
    else:
        # make new working directory, enter it, set up symlinks for the files we need
        hdf5_path = os.path.abspath(fname + '.hdf5') # get /full/path/to/fname.hdf5
        old_working_directory = os.getcwd() # get the current working directory
        os.mkdir(dir_path) # create the output directory
        os.chdir(dir_path) # make it the working directory
        os.symlink(hdf5_path, fname + '.hdf5') # create symlink to .hdf5 file
        os.system('ln -s ~/gizmo-files/* ./') # create symlinks to other GIZMO files that are needed

        # extract data from .hdf5 file
        gas_fraction, disk_mass, memory, memory_total, time_limit_in_s, time_hms = get_mem_time_info(1e7, hdf5_path, n_tasks)

        # create the GIZMO job batch script
        with open(fname + '.gizmojob.sh', 'w') as f:
            f.write(template_slurm_script.format(fname, exec_path, memory_total, time_hms, n_tasks))

        # create the GIZMO parameter file
        with open(fname + '.param', 'w') as f:
            f.write(template_param_file.format(fname, memory, time_limit_in_s, gas_fraction, disk_mass))

        # schedule the job with slurm
        os.system('sbatch ' + fname + '.gizmojob.sh')  

        # go back to the old working directory
        os.chdir(old_working_directory)


