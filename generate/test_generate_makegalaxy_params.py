# this script serves as a test to make sure that generate_makegalaxy_params.py is working
# we do two things:
# 1. run the main function for the parameters Doug used in his 2020 paper
# 2. show a diff of the generated files to the expected output so we can inspect them
#    to make sure they're the same

# PART 1

from generate_makegalaxy_params import generate_makegalaxy_params

gas_mass_resol = 1.e7
output_folder = './makegalaxy_test'
mass_file = './masses_Rennehan_2020.csv'

generate_makegalaxy_params(gas_mass_resol, output_folder, mass_file, REDSHIFT=4.434, Omega_m0=0.3, Omega_L0=0.7, verbose=True)

# PART 2

import os

# folder containing the expected output
comparison_folder = './makegalaxy_Rennehan_Sept2018_1e7'

# display a diff
os.system('diff {} {} | less'.format(output_folder, comparison_folder))

