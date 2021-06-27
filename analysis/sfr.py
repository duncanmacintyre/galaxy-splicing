import pandas
import matplotlib.pyplot as plt

code_mass_units_in_Msun = 1e10
code_time_units_in_Myr = 978.028

df = pandas.read_csv('sfr.txt', sep=' ', names=('t', 'expected', 'SFR', 'M* per step', 'M* this step'))

df['t'] = df['t'] * code_time_units_in_Myr
df['SFR'] = df['SFR'] * code_mass_units_in_Msun / (code_time_units_in_Myr * 1e6)

fig = plt.figure()
ax = plt.axes()
df.plot('t', 'SFR', ax=ax)
ax.set_xlabel('Time [Myr]')
ax.set_ylabel('SFR [M☉/yr]')
fig.savefig('sfr.pdf')
