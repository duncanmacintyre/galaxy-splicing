# DM: wrote June 2021
#
#   plot_mass_envelopes.py
#
#   This file includes a function to plot:
#       * How the stellar mass changes over time;
#       * How the stellar mass envelope (mass within certain radii) changes over time; and
#       * How the radii containing half and 3/4 of the stellar mass change over time.
#

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm
import matplotlib.colors

# plot a figure showing how the a galaxy's mass changes over time and where it's found
#
#   t is a 1D array of time indices
#   total is a 1D array of the total stellar mass (over time)
#   mass_matrix is a 2D array of mass amounts within set radii: cols for different radii, rows for t
#   radii is a 1D array of radii for each of the columns in mass_matrix (determines line colours)
#   t_units is a string of time units
#   mass_units is a string of mass units
#   radii_units is a string of radius units
#
#   title is a string to use for the axes title; use None for no title
#   mass_name is a string for axis labels, e.g. "stellar mass", "gas mass", ...; default just "mass"
#   if relative=True, will plot fraction of mass contained rather than absolute amount
#   if radii_legend=True, will draw legend with the coloured lines for different radii
#   if colorbar=True, will draw colorbar for line colours
#
# returns fig, ax where fig is a Figure object and ax is an array of Axes objects
#
def plot_mass_envelopes(t, total, mass_matrix, radii, t_units, mass_units, radii_units, title=None,
                        mass_name="mass", relative=False, radii_legend=False, colorbar=True):

    total = np.array(total)
    radii = np.array(radii)
    if relative:
        mass_matrix = (np.array(mass_matrix).T/total).T
        # formed = np.array(formed)/total
        total = np.ones(total.shape)
    else:
        mass_matrix = np.array(mass_matrix)
        # formed = np.array(formed)

    # plot of how the stellar mass changes over time, showing the fraction formed since start
    # fig1, ax = plt.subplots()
    # ax.plot(t, formed, color="lightgreen", label="formed")
    # ax.plot(t, total, color="darkblue", label="total")
    # ax.legend()
    # ax.set_ylabel("Stellar mass [{}]".format(mass_units))
    # ax.set_xlabel("Time since start [{}]".format(t_units))

    # get a ColorMap instance
    cm = matplotlib.cm.get_cmap("viridis_r")
    # get a Normalize instance; when called, this will map radii onto the interval [0, 1]
    normalize = matplotlib.colors.Normalize(vmin=0, vmax=np.max(radii))

    # plot of how the stellar mass envelope (mass within certain radii) changes over time
    fig, ax1 = plt.subplots(figsize=(9,6), constrained_layout=True)
    ax2 = ax1.twinx() # this overlaps with ax1 but has a different scale (on the right)
    for y, r in zip(mass_matrix.T, radii):
        ax1.plot(t, y, color=cm(normalize(r)), linewidth=0.5,
                 label="{:.1f} {}".format(r, radii_units))
    if not relative:
        ax1.plot(t, total, color="k", label="total")
    # ax1.plot(t, formed, color="cyan", linestyle="--", label="formed")
    ax2.plot(t, compute_radius_within_frac_mass(total, mass_matrix, radii, frac=0.75),
              color="brown", linestyle=":", label="0.75")
    ax2.plot(t, compute_radius_within_frac_mass(total, mass_matrix, radii, frac=0.5),
              color="red", linestyle=":", label="0.5")
    if relative:
        ax1.set_ylabel("Fraction of {} contained within cylinder".format(mass_name))
    else:
        ax1.set_ylabel("{} contained within cylinder [{}]".format(firstcap(mass_name), mass_units))
    ax1.set_xlabel("Time since start [{}]".format(t_units))
    ax2.set_ylabel("Radius containing fraction of {} [{}]".format(mass_name, radii_units))
    # ensure vertical scales start at 0
    if relative:
        ax1.set_ylim([0, 1])
    else:
        ax1.set_ylim([0, ax1.get_ylim()[1]])
    ax2.set_ylim([0, ax2.get_ylim()[1]])
    # remove extra space from x scale
    ax1.set_xlim([np.min(t), np.max(t)])
    if radii_legend:
        ax1.legend(bbox_to_anchor=(-0.115, 0.5), loc='center right', borderaxespad=0., ncol=3)
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
def compute_radius_within_frac_mass(total, mass_matrix, radii, frac):
    return(np.fromiter(
        (np.interp(targetmass, m, radii) for targetmass, m in zip(frac*total, mass_matrix)),
        dtype="float",
        count=len(total)
    ))

# given a string, return the same but with the first character as uppercase
def firstcap(s):
    return(s[0].upper() + s[1:])


