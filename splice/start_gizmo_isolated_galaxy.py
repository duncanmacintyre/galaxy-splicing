# this script schedules a GIZMO job to simulate a galaxy in isolation
# we do this to ensure structural stability before splicing galaxies together

import os
import sys
import math

import h5py
import numpy as np

from common import grab_property, code_mass_to_Msun


# ----- constants -----

# in this directory, we will create a simulation subdirectory for each galaxy
gizmo_top_dir = os.path.abspath('/scratch/{}/gizmo-galaxies-in-isolation-1e7'.format(os.getenv('USER')))

# how many CPUs are there on each node on this cluster?
cpus_on_node = 48
# how much memory is there on each node on this cluster, in MB?
mem_on_node = 4010*cpus_on_node

# paths to Slurm script template and GIZMO parameter file template
# os.path.dirname(__file__) is the folder containing this script
path_to_slurm_script_template = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                             'protocluster.slurm_job.sh')
path_to_gizmo_param_template = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                            'isolation-1e7.param')

# which version of GIZMO should we use? this is the commit SHA
# we expect the binary to be at ~/bin/GIZMO_<commit_SHA>
gizmo_commit = '9b787c1'

# email address for slurm email notifications
email = ''


# ----- functions -----

# given a .hdf5 output file from makegalaxy, return particle amounts and suggest resources
#
# arguments:
#       * path to .hdf5 file to analyze
#
# returns tuple of ten strings:
#       * gas fraction
#       * mass per halo particle (code units)
#       * total disk mass (in Msun)
#       * total halo mass (in Msun)
#       * suggested Gizmo memory limit per task based on number of particles (in MB)
#       * suggested total memory to request per node in the Slurm script (in MB)
#       * suggested time limit in seconds
#       * suggested time limit in hh:mm:ss form
#       * suggested number of nodes (either 1 or 4)
#       * suggested number of processes per node (either 4, 16 or cpus_on_node)
#
# We assume that (1) there are no bulge or star particles, and (2) we should get masses 
# from the mass table in the header, not from individual particle masses.
#
def suggest_resources(fname):

    with h5py.File(fname, 'r') as f:
        mass_table = np.asarray(f['/Header'].attrs['MassTable']).reshape((-1,))
        N_gas = f['/PartType0/Coordinates'].shape[0]
        N_disk = f['/PartType2/Coordinates'].shape[0]
        N_halo = f['/PartType1/Coordinates'].shape[0]

    # ---- gas fraction ----
    gas_fraction = '{:.8f}'.format(N_gas*mass_table[0]/(N_gas*mass_table[0]+N_disk*mass_table[2]))

    # ---- mass per halo particle (in code units) ----
    halo_res = '{:.9g}'.format(mass_table[1])

    # ---- total disk and halo mass (in Msun) ----
    disk_mass = '{:.9g}'.format(code_mass_to_Msun(N_disk * mass_table[2]))
    halo_mass = '{:.9g}'.format(code_mass_to_Msun(N_halo * mass_table[1]))

    # ---- suggested number of tasks
    x = N_gas + N_disk + N_halo # total number particles
    n_nodes = 4 if x > 400000 else 1
    if x > 30000:
        ntasks_per_node = cpus_on_node
    elif x > 2000:
        ntasks_per_node = 16
    else:
        ntasks_per_node = 4
    n_nodes_string = '{:d}'.format(n_nodes)
    ntasks_per_node_string = '{:d}'.format(ntasks_per_node)

    # ---- suggested memory per task based on number of particles (in MB) ---
    if ntasks_per_node == cpus_on_node: # if we're using whole nodes, ask for all memory
        memory = mem_on_node/cpus_on_node - 1 # memory available per task
        memory_total = 0 # we use 0 to request all the memory on the node
    elif ntasks_per_node < cpus_on_node: # if we're only using part of a node, estimate how much memory needed
        memory = suggest_run_memory(x, n_nodes*ntasks_per_node) # in MB, per task
        memory_total = memory*ntasks_per_node + 1 # in MB, per node
    else: # ntasks_per_node > cpus_on_node; this should never happen!
        assert(False)
    memory_string = '{:.0f}'.format(memory) # per task
    memory_total_string = '{:.0f}'.format(memory_total) # total 

    # ---- suggested time limit ----
    time_limit_in_s = suggest_run_time(x, n_nodes*ntasks_per_node) # in s
    time_string = str(time_limit_in_s) # in s
    # from https://stackoverflow.com/a/775075/13326516
    m, s = divmod(time_limit_in_s, 60)
    h, m = divmod(m, 60)
    time_hms = '{:02d}:{:02d}:{:02d}'.format(h, m, s) # in hh:mm:ss form

    return (gas_fraction, mass_table[1], disk_mass, halo_mass, memory_string, memory_total_string,
            time_string, time_hms, n_nodes_string, ntasks_per_node_string)


# given number of particles x and number of processors p, return a suggested simulation run time in seconds
def suggest_run_time(x, p):
    return int(60*45 + (0.3 * x)/math.sqrt(p))

# given number of particles x and number of processors p, return a suggested per-task memory limit in MB
a = 300 # horizontal scaling factor
m = 2500 # min possible memory value we might assign per task (in MB)
M = 4100 # max possible memory value we might assign per task (in MB)
def suggest_run_memory(x, p):
    return max(
                 ((M-500)/1500)
                 * (
                     math.log10((a*x/(p*(M-500)) + 400)/100000)
                     / math.log10((a*x/(p*(M-500)) + 400)/100000 + 1)
                     + 1499
                   )
                 + 500,


                 m
        )


# ----- beginning of script part -----


# load templates for Gizmo parameter file and Slurm job script
with open(path_to_slurm_script_template) as f:
    template_slurm_script = f.read()
with open(path_to_gizmo_param_template) as f:
    template_param_file = f.read()

# iterate over .hdf5 files; fname is the file name without path or extension
for fname in (os.path.splitext(os.path.basename(s))[0] for s in sys.argv[1:]):

    # path to directory in which we'll store outputs
    dir_path = os.path.join(gizmo_top_dir, fname)

    if os.path.exists(dir_path):
        print('WARNING: {} already exists, skipping {}'.format(dir_path, fname))
    else:
        # make new working directory, enter it, set up symlinks for the files we need
        hdf5_path = os.path.abspath(fname + '.hdf5') # get /full/path/to/fname.hdf5
        old_working_directory = os.getcwd() # get the current working directory
        os.mkdir(dir_path) # create the output directory
        os.chdir(dir_path) # make it the working directory
        os.symlink(hdf5_path, fname + '.hdf5') # create symlink to .hdf5 file
        os.system('ln -s ~/gizmo-files/* ./') # create symlinks to other GIZMO files that are needed

        # extract data from .hdf5 file
        (gas_fraction, halo_res, disk_mass, halo_mass, memory, memory_total,
            time_limit_in_s, time_hms, n_nodes, ntasks_per_node) = suggest_resources(hdf5_path)

        formatter = {
                'ic_file_name':fname, # we want no extension for this
                'param_file_name':fname + '.param',
                'job_name':fname,
                'resubmit_slurm_script':fname + '.gizmojob.sh',
                'timelim_string':time_hms,
                'timelim_in_s':time_limit_in_s,
                'nodes':n_nodes,
                'ntasks_per_node':ntasks_per_node,
                'memory_per_node':memory_total,
                'memory_per_task':memory,
                'resubmit_on':'1',
                'sim_time_max':'0.2',
                'black_hole_res':halo_res,
                'gas_fraction':gas_fraction,
                'stellar_mass':disk_mass,
                'halo_mass':halo_mass,
                'email':email,
                'gizmo_commit':gizmo_commit
            }

        # create the GIZMO job batch script
        with open(fname + '.gizmojob.sh', 'w') as f:
            f.write(template_slurm_script.format_map(formatter))

        # create the GIZMO parameter file
        with open(fname + '.param', 'w') as f:
            f.write(template_param_file.format_map(formatter))

        # schedule the job with slurm
        os.system('sbatch ' + fname + '.gizmojob.sh')  

        # go back to the old working directory
        os.chdir(old_working_directory)


