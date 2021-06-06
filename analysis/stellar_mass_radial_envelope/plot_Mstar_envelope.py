#!/bin/python
#
#   plot_Mstar_envelope.py
#
#   The file org_fgas=0.5_star_mass_data.txt contains stellar mass data from Doug's original
#   simulation with a gas fraction of 0.5. This script produces three plots of this data:
#       * How the stellar mass changes over time, showing the fraction formed since start;
#       * How the stellar mass envelope (mass within certain radii) changes over time; and
#       * How r_(1/2), the radius containing half the stellar mass, changes over time.
#

import numpy as np
import pandas
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap



# plot a figure showing how the stellar mass changes over time and where it's found
#
#   t is a 1D array of time indices
#   total is a 1D array of the total stellar mass (over time)
#  *formed is a 1D array of the stellar mass formed since the start of the simulation (over time)
#   mass_matrix is a 2D array with columns each giving the stellar mass within a certain radius and
#       and rows corresponding to rows for the time indices
#   radii is a 1D array of radii for each of the columns in mass_matrix (determines line colours)
#   radii_names is a 1D array of strings giving names for the radii (to include in legends)
#   t_units is a string of time units
#   mass_units is a string of stellar mass units
#   radii_units is a string of radius units
#
#   * not currently in use
#
# returns fig, ax where fig is a Figure object and ax is an array of Axis objects
#
def plotMassEnvelopes(t, total, formed, mass_matrix, radii, radii_names,
                      t_units, mass_units, radii_units, relative=False):

    # get a ColorMap instance
    cm = get_cmap()

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
    fig, ax1 = plt.subplots(figsize=(9,5), constrained_layout=True)
    ax2 = ax1.twinx() # this overlaps with ax1 but has a different scale (on the right)
    colors = [cm(x) for x in 1-radii/np.max(radii)]
    for y, c, s in zip(mass_matrix.T, colors, radii_names):
        ax1.plot(t, y, color=c, label=s)
    if not relative:
        ax1.plot(t, total, color="k", label="total")
    # ax1.plot(t, formed, color="cyan", linestyle="--", label="formed")
    ax2.plot(t, computeRadiusWithFracMass(total, mass_matrix, radii, frac=0.75),
              color="brown", linestyle=":", label="0.75")
    ax2.plot(t, computeRadiusWithFracMass(total, mass_matrix, radii, frac=0.5),
              color="red", linestyle=":", label="0.5")
    if relative:
        ax1.set_ylabel("Fraction of stellar mass within radius [{}]".format(mass_units))
    else:
        ax1.set_ylabel("Stellar mass within radius [{}]".format(mass_units))
    ax1.set_xlabel("Time since start [{}]".format(t_units))
    ax2.set_ylabel("Radius containing fraction of mass [{}]".format(radii_units))
    # ensure vertical scales start at 0
    ax1.set_ylim([0, ax1.get_ylim()[1]])
    ax2.set_ylim([0, ax2.get_ylim()[1]])
    # remove extra space from x scale
    ax1.set_xlim([np.min(t), np.max(t)])
    ax1.legend(bbox_to_anchor=(-0.115, 0.5), loc='center right', borderaxespad=0.)
    ax2.legend(bbox_to_anchor=(1.115, 0.5), loc='center left', borderaxespad=0.)
    
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
data = pandas.read_csv("org_fgas=0.5_star_mass_data.txt", sep="\t")
# radii at which data binned - must be increasing order
radii = np.array([5., 10., 15., 30., 50., 70.])
# columnn names for radial bins in data
radii_names = ["{:.0f}kpc".format(n) for n in radii]

# make the plots
fig1, _ = plotMassEnvelopes(data["Time"], data["Total"], data["Formed"], data[radii_names],
                       radii, radii_names, "Gyr", "M☉", "kpc", relative=False)
fig2, _ = plotMassEnvelopes(data["Time"], data["Total"], data["Formed"], data[radii_names],
                       radii, radii_names, "Gyr", "M☉", "kpc", relative=True)

fig1.savefig("org_fgas=0.5_absolute_Mstar_envelope.pdf")
fig2.savefig("org_fgas=0.5_relative_Mstar_envelope.pdf")

plt.show()

