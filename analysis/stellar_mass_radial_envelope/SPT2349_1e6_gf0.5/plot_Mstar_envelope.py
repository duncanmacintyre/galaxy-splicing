import sys, os, pandas
import numpy as np
import matplotlib.pyplot as plt
# add parent directory to path, per https://stackoverflow.com/a/30536516/13326516
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# import plotting function from parent directory
from plot_mass_envelopes import plot_mass_envelopes

def plot_and_save(name, file, prefix):

    # load the data
    data = pandas.read_fwf(file)
    # radii at which data binned - must be increasing order
    radii = np.loadtxt("../../radius_bins.txt")
    # radii = np.array([5., 10., 15., 30., 50., 70.])
    radii_names = ["{}{:.1f}kpc".format(prefix, n) for n in radii]

    # make the plots
    fig1, _ = plot_mass_envelopes(data["Time"], data["Total"], data[radii_names], radii, 
                                  "Myr", "M☉", "kpc", title="$f_{gas} = 0.5$", mass_name="stellar mass", relative=False)
    fig2, _ = plot_mass_envelopes(data["Time"], data["Total"], data[radii_names], radii,
                                  "Myr", "M☉", "kpc", title="$f_{gas} = 0.5$", mass_name="stellar mass", relative=True)

    fig1.savefig("org_fgas=0.5_{}_absolute_Mstar_envelope.pdf".format(name))
    fig2.savefig("org_fgas=0.5_{}_relative_Mstar_envelope.pdf".format(name))

plot_and_save("mean", "star_mass_data_mean.txt", "")
plot_and_save("x", "star_mass_data.txt", "x")
plot_and_save("y", "star_mass_data.txt", "y")
plot_and_save("z", "star_mass_data.txt", "z")

plt.show()