
from generate_makegalaxy_params import generate_makegalaxy_params

generate_makegalaxy_params(
        gas_mass_resol=1.e6,
        output_folder='./makegalaxy-1e6',
        mass_file='./galaxy_masses_105.txt',
        inflate_fgas={'C9': 0.25}
    )

generate_makegalaxy_params(
        gas_mass_resol=1.e7,
        output_folder='./makegalaxy-1e7',
        mass_file='./galaxy_masses_105.txt',
        inflate_fgas={'C9': 0.25}
    )
