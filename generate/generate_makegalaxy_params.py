# DR 2020: intial code for 2020 paper
# DM 2021: reorganized code; included more galaxies, new Mstars data
#
# The file defines the function generate_makegalaxy_params(). This function
# generates parameter files to use as input to makegalaxy based on galaxy
# gas masses, stellar masses, and halo virial masses.
#
# I recommend importing generate_makegalaxy_params into a script and running
# it from there. That way you keep track of all the options that you used.
# See for example gen_files.py.
#
# The python script 'test_generate_makegalaxy_params.py' can be run to test
# that the function is working.
#

import numpy as np
import pandas
import itertools

### CONFIGUATION OPTIONS ###

# conversion from data file units to mass of sun
fac = 1.0e10

# How much bigger should the DM/Bulge particles be than gas/stars?
resol_factor = 5.0

h = 0.7
G = 4.302e-9 # Mpc (km/s)**2 M_sun**-1
H0 = 100 # h km Mpc**-1 s**-1    # Hubble's constant v = H0 D
M_collapse = 8e12 # Msun h**-1

keys = ['OutputDir',
        'OutputFile',
        'CC',
        'V200',
        'LAMBDA',
        'AnisotropyRadius',
        'MD',
        'MB',
        'MBH',
        'N_HALO',
        'N_DISK',
        'N_GAS',
        'N_BULGE',
        'UseQrng',
        'WriteDensity',
        'JD',
        'H',
        'DiskHeight',
        'RadialDispersionFactor',
        'DiskPopMode',
        'DiskPopAge',
        'DiskPopTau',
        'GasFraction',
        'MaxGasDiskHeight',
        'GasDistribution',
        'GasExpAlpha',
        'PowerLawGamma',
        'PowerLawCutoff',
        'Zmet0',
        'ZmetGrad',
        'GravSoftening',
        'ErrTolTheta',
        'BulgeSize',
        'BulgePopMode',
        'BulgePopAge',
        'BulgePopTau',
        'MaxSfrTimescale',
        'FactorSN',
        'FactorEVP',
        'TempSupernova',
        'TempClouds',
        'FactorForSofterEQS',
        'REDSHIFT',
        'Omega_m0',
        'Omega_L0']

# default values for the galaxy, some will be overwritten
values = ['./',
            '',
            '3.5',                      # Halo Concentration
            '230.0',                    # V200
            '0.033',                    # Spin parameter lambda
            '0.0',                      # Anisotropy radius
            '0.064',                    # Disk mass fraction
            '0.045',                    # Bulge mass fraction
            '0.000000038',              # Black hole mass
            '1750000',                  # Num. of halo particles
            '375000',                   # Num. of disk particles
            '875000',                   # Num. of gas particles
            '875000',                   # Num. of bulge particles
            '0',                        # 0=pseudorandom, 1=quasirandom
            '0',                        # write density instead of mass
            '0.064',                    # Disk spin fraction
            '0',                        # Disk scale L, usually set by disk spin
            '0.1',                      # SD height in units of radial SL
            '1.0',                      # Radial dispersion
            'exponential',              # Stellar population mode
            '13.9',                     # Disk population age (in Gyr) - doesn't matter, age has no effect
            '-106.0',                   # Disk population tau (in Gyr)
            '0.7',                      # Gas fraction
            '2.0',                      # Max disk height
            '0',                        # Gas distribution type
            '1.0',                      # GD scale length multiplier
            '1.0',                      # Power law gamma
            '1.0',                      # Max. radius of gas for PLD
            '0.032',                    # Gas metallicity in solar (central)
            '-0.03',                    # Gas metallicity gradient (dex/kpc)
            '0.06',                     # Grav. softening
            '0.15',                     # Tree opening criterion
            '0.75',                     # Bulge SL in units of disk SL
            'instantaneous',            # Stellar population mode (bulge)
            '13.9',                     # Bulge population age (in Gyr)
            '1.0',                      # Bulge population tau (in Gyr)
            '4.5',                      # Maxsfrtimescale
            '0.1',                      # Beta, mass fraction of massive stars
            '3000',                     # A_0 for SH+03 model
            '3.0e8',                    # T_SN effective SN temperature
            '1000',                     # Temperature of cold clouds
            '1.0',                      # q_eos
            '4.303',                    # redshift
            '0.3089',                   # omega_matter
            '0.6911']                   # omega_lambda

# build mapping for keys and values
idx = {}
for i, key in enumerate(keys):
    idx.update({key: i})


### MAIN FUNCTION TO GENERATE PARAMTER FILES ###

# generate parameter files for makegalaxy
#       gas_mass_resol is the mass per gas particle (units: mass of the sun)
#       output_folder is the path in which to save the parameter files
#       mass_file is the path to a CSV file with columns label, Mgas, M*, Mvir
#       write is whether to actually write any files (default True, of course)
#       verbose is whether to print particle numbers for each galaxy (default False)
# You can pass additional keyword arguments. These will become new default values for
# makegalaxy parameters. e.g. pass REDSHIFT=2 to set the redshift parameter to 2.
def generate_makegalaxy_params(gas_mass_resol: float, output_folder: str, mass_file: str,
                               write: bool = True, verbose: bool = False,
                               **kwargs):

    # mass per particle of each kind
    stars_mass_resol = gas_mass_resol
    dm_mass_resol = resol_factor * gas_mass_resol
    bulge_mass_resol = resol_factor * gas_mass_resol

    # --- load galaxy masses from file ---

    # load the gas, stellar, and halo masses from a CSV file
    # the CSV must use commas for delimeters and have at least these columns:
    #   label, Mgas, M*, Mvir
    masses = pandas.read_csv(mass_file) # type: pandas.DataFrame
    # masses['label'] has the labels C1, C2, C3, ... as strings
    # masses['Mgas'] has the gas masses
    # masses['M*'] has the stellar masses
    # masses['Mvir'] has the virulent masses of the dark matter halos

    # which galaxies in masses we will generate parameter files for
    # keep only galaxies where Mgas, M*, and Mvir are not NaN
    # this ignores any galaxies for which we did not compute mass
    keep = masses['Mgas'].notna() & masses['M*'].notna() & masses['Mvir'].notna()

    # number of parameter files that will be generated (one for each galaxy)
    M = keep.sum()

    # this iterator continuously cycles over the galaxy labels
    labels = itertools.cycle(masses.loc[keep, 'label'])
    
    # these masses in units of the mass of the sun
    gas_mass = masses.loc[keep, 'Mgas'].to_numpy() * fac
    stars_mass = masses.loc[keep, 'M*'].to_numpy() * fac
    dm_mass = masses.loc[keep, 'Mvir'].to_numpy() * fac

    # --- compute derived quantities ---

    tot_mass = gas_mass + stars_mass + dm_mass

    concentrations = halo_concentration(tot_mass)
    circ_vels = mass_to_vel(tot_mass)

    disk_mass = gas_mass + stars_mass 
    bulge_mass = np.zeros(M) # assume the bulge could not have formed in 1.4Gyr without many mergers
    disk_mass_frac = disk_mass / tot_mass
    bulge_mass_frac = bulge_mass / tot_mass

    bh_mass_frac = gas_mass_resol / tot_mass # black hole mass is same as mass of one gas particle

    # compute number of particles
    Ngas = gas_mass / gas_mass_resol
    Ndm = dm_mass / dm_mass_resol
    Nstars = stars_mass / stars_mass_resol
    Nbulge = bulge_mass / bulge_mass_resol

    # --- write the parameter files ---

    if write:
        # loop through each row from the data and make the new parameters
        for i in range(0, M):
            if verbose:
                print('Preparing to write for {}'.format(masses.loc[keep, 'label'][i]))
                print('Gas mass, #: %g, %d' % (gas_mass[i], Ngas[i]))
                print('Disk mass, #: %g, %d' % (disk_mass[i], Nstars[i]))
                print('Bulge mass, #: %g, %d' % (bulge_mass[i], Nbulge[i]))
                print('DM mass, #: %g, %d' % (dm_mass[i], Ndm[i]))

            tmp_values = values
            tmp_label = next(labels)

            # overwrite any values passed
            for key_given, value_given in kwargs.items():
                tmp_values[idx[key_given]] = str(value_given)
            
            tmp_values[idx['OutputFile']] = tmp_label + '.hdf5'
            
            tmp_values[idx['CC']] = str(concentrations[i])
            tmp_values[idx['V200']] = str(circ_vels[i])
            
            tmp_values[idx['MB']] = str(float(bulge_mass_frac[i]))
            tmp_values[idx['MD']] = str(float(disk_mass_frac[i]))
            tmp_values[idx['JD']] = str(float(disk_mass_frac[i]))

            tmp_values[idx['N_HALO']] = str(int(Ndm[i]))
            tmp_values[idx['N_DISK']] = str(int(Nstars[i]))
            tmp_values[idx['N_GAS']] = str(int(Ngas[i]))
            tmp_values[idx['N_BULGE']] = str(int(Nbulge[i]))
           
            tmp_values[idx['MBH']] = str(float(bh_mass_frac[i]))

            if verbose:
                print('Writing file ' + output_folder + '/' + tmp_label + '.txt')
                print()
            with open(output_folder + '/' + tmp_label + '.txt', 'w') as f:
                for key in keys:
                    f.write(key + '\t\t\t\t' + tmp_values[idx[key]] + '\n')
        if verbose:
            print('Done writing files')
            print()
    
    # --- compute total number of particles in all galaxies, print this ---

    total_mass = np.sum(tot_mass)

    total_Ngas = np.sum(Ngas)
    total_Ndm = np.sum(Ndm)
    total_Nstars = np.sum(Nstars)
    total_Nbulge = np.sum(Nbulge)

    N = Ngas + Ndm + Nstars + Nbulge
    total_part = np.sum(N)

    if verbose:
        print()
        print('SUMMARY OF PARTICLE NUMBERS')
        print()
        print(pandas.DataFrame({'label':masses.loc[keep, 'label'],
                                'mass':tot_mass,
                                'Ngas':Ngas,
                                'Ndm':Ndm,
                                'Nstars':Nstars,
                                'Nbulge':Nbulge,
                                'total':N}).to_string())
        print()

    print('{} galaxies of total mass {:.3e}'.format(M, total_mass))
    print('TotNgas: %g' % total_Ngas)
    print('TotNdm: %g' % total_Ndm)
    print('TotNstars: %g' % total_Nstars)
    print('TotNbulge: %g' % total_Nbulge)
    print('TotPart: %g' % total_part)

    return(total_mass, total_part, total_Ngas, total_Ndm, total_Nstars, total_Nbulge)


### HELPER FUNCTIONS TO COMPUTE QUANTITIES ###
# we compute all quantities as if the galaxies are at redshift 0
# makegalaxy will then convert quantities to be correct for our non-zero redshift
#
# note:
#       we want to compile makegalaxy with the REDSHIFT_SCALING option set
#       and the V_SCALING option unset in "config.h"
#       (this is the default at time of writing)

# calculate the halo concentration (from Roberston+06)
def halo_concentration(mass):
    return 9.0 * (mass * h / M_collapse)**(-0.13)

# calculate the virial velocity V200
# mass must be in Msun
def mass_to_vel(mass):
    # (mass / h) to follow Springel+05 Section 2.4
    return (10.0 * G * h * H0 * mass)**(1.0 / 3.0)


