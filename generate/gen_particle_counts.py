import numpy as np
from generate_makegalaxy_params import generate_makegalaxy_params
import pandas

output_folder = './makegalaxy'

# set format to print in
pandas.set_option('display.float_format', '{:.2e}'.format)

df1 = pandas.DataFrame({
        'gas_mass_resol': np.repeat(np.array([5.e5, 1.e6, 4.e6, 1.e7]), 4),
        'mass_file': np.repeat(np.array(['all', 'C', '210', 'notSPIREc'], ndmin=2), 4, axis=0).flatten()
    })
df2 = pandas.DataFrame((generate_makegalaxy_params(gas_mass_resol, output_folder, 'galaxy_masses_{}.csv'.format(mass_file), write=False)
                        for gas_mass_resol, mass_file in zip(df1['gas_mass_resol'], df1['mass_file'])), 
                       columns=['mass', 'particles', 'Ngas', 'Ndm', 'Nstars', 'Nbulge'])

df = df2[['mass', 'particles']]
df.index = pandas.MultiIndex.from_frame(df1)

df.to_latex('table_particle_no/table_generated.tex')

print(df)