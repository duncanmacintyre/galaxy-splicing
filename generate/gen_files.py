
from generate_makegalaxy_params import generate_makegalaxy_params

gas_mass_resol = 1.e6
output_folder = './makegalaxy'
mass_file = './galaxy_masses_105.txt'
inflate_fgas = {'C9': 0.25}

generate_makegalaxy_params(gas_mass_resol, output_folder, mass_file, inflate_fgas=inflate_fgas)

