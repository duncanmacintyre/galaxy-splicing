import sys, os, pandas
import numpy as np
import matplotlib.pyplot as plt
# add parent directory to path, per https://stackoverflow.com/a/30536516/13326516
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# import plotting function from parent directory
from plot_mass_envelopes import plot_mass_envelopes


# load the data
data_mean = pandas.read_fwf("star_mass_data_mean.txt")
data_x = pandas.read_fwf("star_mass_data_x.txt")
# radii at which data binned - must be increasing order
radii = np.loadtxt("../../radius_bins.txt")
# radii = np.array([5., 10., 15., 30., 50., 70.])
radii_names = ["{:.1f}kpc".format(n) for n in radii]

# make the plots
fig1, _ = plot_mass_envelopes(data_mean["Time"], data_mean["Total"], data_mean[radii_names], radii, 
                              "Myr", "M☉", "kpc", title="$f_{gas} = 0.9$", mass_name="stellar mass", relative=False)
fig2, _ = plot_mass_envelopes(data_mean["Time"], data_mean["Total"], data_mean[radii_names], radii,
                              "Myr", "M☉", "kpc", title="$f_{gas} = 0.9$", mass_name="stellar mass", relative=True)
fig3, _ = plot_mass_envelopes(data_x["Time"], data_x["Total"], data_x[radii_names], radii, 
                              "Myr", "M☉", "kpc", title="$f_{gas} = 0.9$", mass_name="stellar mass", relative=False)
fig4, _ = plot_mass_envelopes(data_x["Time"], data_x["Total"], data_x[radii_names],radii,
                              "Myr", "M☉", "kpc", title="$f_{gas} = 0.9$", mass_name="stellar mass", relative=True)

fig1.savefig("org_fgas=0.9_mean_absolute_Mstar_envelope.pdf")
fig2.savefig("org_fgas=0.9_mean_relative_Mstar_envelope.pdf")
fig3.savefig("org_fgas=0.9_x_absolute_Mstar_envelope.pdf")
fig4.savefig("org_fgas=0.9_x_relative_Mstar_envelope.pdf")

plt.show()

