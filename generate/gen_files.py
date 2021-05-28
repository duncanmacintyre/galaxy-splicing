
from generate_makegalaxy_params import generate_makegalaxy_params

gas_mass_resol = 1.e7
output_folder = './makegalaxy'
mass_file = './galaxy_masses_105.txt'

generate_makegalaxy_params(gas_mass_resol, output_folder, mass_file)

