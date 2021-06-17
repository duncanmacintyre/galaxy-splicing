#!/bin/python
#
#   plot_Mstar_envelope.py
#
#   The file star_mass_data_mean.txt contains stellar mass data from Doug's original
#   simulation with a gas fraction of 0.7. This script produces three plots of this data:
#       * How the stellar mass changes over time, showing the fraction formed since start;
#       * How the stellar mass envelope (mass within certain radii) changes over time; and
#       * How r_(1/2), the radius containing half the stellar mass, changes over time.
#

import numpy as np
import pandas
import matplotlib.pyplot as plt
import matplotlib.cm
import matplotlib.colors



# plot a figure showing how the stellar mass changes over time and where it's found
#
#   t is a 1D array of time indices
#   total is a 1D array of the total stellar mass (over time)
#   mass_matrix is a 2D array with columns each giving the stellar mass within a certain radius and
#       and rows corresponding to rows for the time indices
#   radii is a 1D array of radii for each of the columns in mass_matrix (determines line colours)
#   radii_names is a 1D array of strings giving names for the radii (to include in legends)
#   t_units is a string of time units
#   mass_units is a string of stellar mass units
#   radii_units is a string of radius units
#
#   title is a string to use for the axes titles
#   if relative=True, will plot fraction of stellar mass contained rather than absolute amount
#
# returns fig, ax where fig is a Figure object and ax is an array of Axes objects
#
def plotMassEnvelopes(t, total, mass_matrix, radii, radii_names,
                      t_units, mass_units, radii_units, title=None, relative=False):

    # get a ColorMap instance
    cm = matplotlib.cm.get_cmap("viridis_r")
    # get a Normalize instance; when called, this will map radii onto the interval [0, 1]
    normalize = matplotlib.colors.Normalize(vmin=np.min(radii), vmax=np.max(radii))

    # colors = [cm(norm(x)) for x in radii]

    # plot of how the stellar mass changes over time, showing the fraction formed since start
    # fig1, ax = plt.subplots()
    # ax.plot(t, formed, color="lightgreen", label="formed")
    # ax.plot(t, total, color="darkblue", label="total")
    # ax.legend()
    # ax.set_ylabel("Stellar mass [{}]".format(mass_units))
    # ax.set_xlabel("Time since start [{}]".format(t_units))

    total = np.array(total)
    radii = np.array(radii)
    if relative:
        mass_matrix = (np.array(mass_matrix).T/total).T
        # formed = np.array(formed)/total
        total = np.ones(total.shape)
    else:
        mass_matrix = np.array(mass_matrix)
        # formed = np.array(formed)

    # plot of how the stellar mass envelope (mass within certain radii) changes over time
    fig, ax1 = plt.subplots(figsize=(9,6), constrained_layout=True)
    # fig = plt.figure(figsize=(9,5), constrained_layout=True)
    # ax1 = fig.add_subplot(5, 1, (1,3))
    ax2 = ax1.twinx() # this overlaps with ax1 but has a different scale (on the right)
    # colors = [cm(x) for x in 1-radii/np.max(radii)]
    for y, r, s in zip(mass_matrix.T, radii, radii_names):
        ax1.plot(t, y, color=cm(normalize(r)), label=s, linewidth=0.5)
    if not relative:
        ax1.plot(t, total, color="k", label="total")
    # ax1.plot(t, formed, color="cyan", linestyle="--", label="formed")
    ax2.plot(t, computeRadiusWithFracMass(total, mass_matrix, radii, frac=0.75),
              color="brown", linestyle=":", label="0.75")
    ax2.plot(t, computeRadiusWithFracMass(total, mass_matrix, radii, frac=0.5),
              color="red", linestyle=":", label="0.5")
    if relative:
        ax1.set_ylabel("Fraction of stellar mass contained within cylinder [{}]".format(mass_units))
    else:
        ax1.set_ylabel("Stellar mass contained within cylinder [{}]".format(mass_units))
    ax1.set_xlabel("Time since start [{}]".format(t_units))
    ax2.set_ylabel("Radius containing fraction of mass [{}]".format(radii_units))
    # ensure vertical scales start at 0
    if relative:
        ax1.set_ylim([0, 1])
    else:
        ax1.set_ylim([0, ax1.get_ylim()[1]])
    ax2.set_ylim([0, ax2.get_ylim()[1]])
    # remove extra space from x scale
    ax1.set_xlim([np.min(t), np.max(t)])
    # ax1.legend(bbox_to_anchor=(0.115, -0.9), loc='center right', borderaxespad=0., ncol=3)
    ax2.legend(bbox_to_anchor=(1.115, 0.5), loc='center left', borderaxespad=0.)
    
    if title is not None:
        ax1.set_title(title)

    fig.colorbar(
        matplotlib.cm.ScalarMappable(norm=normalize, cmap=cm),
        ax=[ax1, ax2], location="bottom", aspect=35, label="Radius [{}]".format(radii_units))

    return(fig, [ax1, ax2])

# for some object, compute the radius containing a fraction of the mass using linear interpolation
#   total is a 1D array of the total mass (e.g., over time, etc.)
#   mass_matrix is a 2D array with columns each giving the stellar mass within certain radius and
#       and rows corresponding to rows in total
#   radii is a 1D array of radii for each of the columns in mass_matrix
#   frac is a float giving the fraction for which we want to compute the radii for; default 0.5
# returns a 1D array the same shape as total
def computeRadiusWithFracMass(total, mass_matrix, radii, frac=0.5):
    return(np.fromiter(
        (np.interp(targetmass, m, radii) for targetmass, m in zip(frac*total, mass_matrix)),
        dtype="float",
        count=len(total)
    ))


# ----- beginning of script -----

# load the data
data_mean = pandas.read_fwf("star_mass_data_mean.txt")
data_x = pandas.read_fwf("star_mass_data_x.txt")
# radii at which data binned - must be increasing order
radii = np.loadtxt("../../radius_bins.txt")
# radii = np.array([5., 10., 15., 30., 50., 70.])
# columnn names for radial bins in data
radii_names = ["{:.1f}kpc".format(n) for n in radii]

# make the plots
fig1, _ = plotMassEnvelopes(data_mean["Time"], data_mean["Total"], data_mean[radii_names], radii, radii_names, 
                       "Myr", "M☉", "kpc", title="$f_{gas} = 0.7$ real01", relative=False)
fig2, _ = plotMassEnvelopes(data_mean["Time"], data_mean["Total"], data_mean[radii_names],radii, radii_names, 
                       "Myr", "M☉", "kpc", title="$f_{gas} = 0.7$ real01", relative=True)
fig3, _ = plotMassEnvelopes(data_x["Time"], data_x["Total"], data_x[radii_names], radii, radii_names, 
                       "Myr", "M☉", "kpc", title="$f_{gas} = 0.7$ real01", relative=False)
fig4, _ = plotMassEnvelopes(data_x["Time"], data_x["Total"], data_x[radii_names],radii, radii_names, 
                       "Myr", "M☉", "kpc", title="$f_{gas} = 0.7$ real01", relative=True)

fig1.savefig("org_fgas=0.7_real01_mean_absolute_Mstar_envelope.pdf")
fig2.savefig("org_fgas=0.7_real01_mean_relative_Mstar_envelope.pdf")
fig3.savefig("org_fgas=0.7_real01_x_absolute_Mstar_envelope.pdf")
fig4.savefig("org_fgas=0.7_real01_x_relative_Mstar_envelope.pdf")

plt.show()

