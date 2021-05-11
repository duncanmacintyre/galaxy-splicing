# DR 2020: intial code for 2020 paper
# DM 2021: modified to include more galaxies, new Msun data

# this file stores the main function
# the python script 'gen_files.py' should be used to run this
# the python script 'test_generate_makegalaxy_params.py' can be run to test that the function is working

# generate parameter files for makegalaxy
#       gas_mass_resol is the mass per gas particle (units: mass of the sun)
#       output_folder is the path in which to save the parameter files
#       mass_file is the path to a CSV file with columns label, Mgas, M*, Mvir
def generate_makegalaxy_params(gas_mass_resol: float, output_folder: str, mass_file: str):

    import numpy as np
    import pandas
    import itertools

    # mass per particle of each kind
    resol_factor = 5.0  # How much bigger should the DM/Bulge particles be?
    stars_mass_resol = gas_mass_resol
    dm_mass_resol = resol_factor * gas_mass_resol
    bulge_mass_resol = resol_factor * gas_mass_resol

    # conversion from data file units to mass of sun
    fac = 1.0e10

    h = 0.7
    G = 4.302e-9 # Mpc (km/s)**2 M_sun**-1
    H0 = 100 # h km Mpc**-1 s**-1    # Hubble's constant v = H0 D
    M_collapse = 8e12 # Msun h**-1

    omega_matter = 0.3
    omega_lambda = 0.7
    # redshift = 4.434
    redshift = 0

    # load the gas, stellar, and halo masses
    masses = pandas.read_csv(mass_file) # type: pandas.DataFrame
    # masses['label'] has the labels C1, C2, C3, ... as strings
    # masses['Mgas'] has the gas masses
    # masses['M*'] has the stellar masses
    # masses['Mvir'] has the virulent masses of the dark matter halos

    # this iterator continuously cycles over the galaxy labels
    labels = itertools.cycle(masses['label'])

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

    # Default values for the galaxy, some will be overwritten
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
                '0',                       # Disk scale L, usually set by disk spin
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
                str(redshift),                   # redshift
                '0.3',                      # omega_matter
                '0.7']                      # omega_lambda

    # Build mapping
    idx = {}
    for i, key in enumerate(keys):
        idx.update({key: i})
        
    def hubble_parameter():
        return H0 * h * np.sqrt(omega_matter * (1.0 + redshift)**3  + omega_lambda)

    # calculate the concentration expected at a redshift
    # from Roberston+06
    def halo_concentration(mass):
        #return 9.0 * (mass / (M_collapse * h))**(-0.13) / (1.0 + redshift)
        return 9.0 * (mass * h / M_collapse)**(-0.13) / (1.0 + redshift)

    # -- this function not used --
    # vel must be in km/s
    # def vel_to_mass(vel):    
    #     return vel**3 / (10.0 * G * hubble_parameter())

    # Mass must be in Msun
    def mass_to_vel(mass):
        # (mass / h) to follow Springel+05 Section 2.4
        #return (10.0 * G * hubble_parameter() * (mass / h))**(1.0 / 3.0)
        return (10.0 * G * hubble_parameter() * mass)**(1.0 / 3.0)

    # these masses in units of the mass of the sun
    gas_mass = masses['Mgas'] * fac
    stars_mass = masses['M*'] * fac
    dm_mass = masses['Mvir'] * fac
    tot_mass = gas_mass + stars_mass + dm_mass

    concentrations = halo_concentration(tot_mass)
    circ_vels = mass_to_vel(tot_mass)

    disk_mass = gas_mass + stars_mass 
    bulge_mass = np.zeros(len(masses)) # assume the bulge could not have formed in 1.4Gyr without many mergers
    disk_mass_frac = disk_mass / tot_mass
    bulge_mass_frac = bulge_mass / tot_mass

    bh_mass_frac = gas_mass_resol / tot_mass # black hole mass is same as mass of one gas particle

    # compute number of particles
    Ngas = gas_mass / gas_mass_resol
    Ndm = dm_mass / dm_mass_resol
    Nstars = stars_mass / stars_mass_resol
    Nbulge = bulge_mass / bulge_mass_resol

    # loop through each row from the data and make the new parameters
    for i in range(0, len(masses)):
        # print('Gas mass, #: %g, %d' % (gas_mass[i], Ngas[i]))
        # print('Disk mass, #: %g, %d' % (disk_mass[i], Nstars[i]))
        # print('Bulge mass, #: %g, %d' % (bulge_mass[i], Nbulge[i]))
        # print('DM mass, #: %g, %d' % (dm_mass[i], Ndm[i]))
        # print()

        tmp_values = values
        tmp_label = next(labels)
        
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

        with open(output_folder + '/' + tmp_label + '.txt', 'w') as f:
            for key in keys:
                f.write(key + '\t\t\t\t' + tmp_values[idx[key]] + '\n')
        

    total_Ngas = np.sum(Ngas)
    total_Ndm = np.sum(Ndm)
    total_Nstars = np.sum(Nstars)
    total_Nbulge = np.sum(Nbulge)

    total_part = total_Ngas + total_Ndm + total_Nstars + total_Nbulge

    print('TotNgas: %g' % total_Ngas)
    print('TotNdm: %g' % total_Ndm)
    print('TotNstars: %g' % total_Nstars)
    print('TotNbulge: %g' % total_Nbulge)
    print('TotPart: %g' % total_part)

