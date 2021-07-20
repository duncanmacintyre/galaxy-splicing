# common functions used in multiple other files
# DM: initial version July 2021

import copy

import numpy as np
from matplotlib.cm import get_cmap
from matplotlib.colors import LogNorm


# --- code_time_to_Myr
# given time in code units, return time in Myr
def code_time_to_Myr(t):
    return 978.028*t


# --- locate_peak_density_3D
# given an Nx3 array a of coordinates, return coordinates of the location with greatest density
#   Bins coordinates into nbins^3 bins and finds bin with greatest count or total weighting.
#   Only considers the domain (-cube_radius, cube_radius) in each of the three directions.
#   Weights can be None or an array of shape (N,).
#
# if return_histogram is False, returns:
#   coordinates_of_maximum        numpy array with coordinates for peak density, shape (3,)
#                                   or numpy.array([0, 0, 0]) if a is empty
#
# if return_histogram is True, returns a length-four tuple with:
#   coordinates_of_maximum        (as above)
#   edges                         numpy array of bin edges, shape (nbins+1,)
#   d_edges                       float with width of bins (they are evenly spaced)
#   unique                        indices for right edges of non-empty bins, shape (M, 3)
#                                   e.g. edges[unique[i,:]] gives the right edges of the ith
#`                                  non-empty bin in the x, y, and z directions
#   counts                        count or total weight inside each non-empty bin, shape (M, 3)
#
def locate_peak_density_3D(a, cube_radius, nbins, weights=None, return_histogram=False):
    # construct array with edges of bins (same edges for each of the three dimensions)
    # there will be nbins^3 bins; each bin can be described by the indices for the right edge:
    #   for example, a bin with indices 100, 20, and 31 could be said to occupy the cube with
    #   edges[99] ≤ x < edges[100], edges[19] ≤ y < edges[20]], and edges[30] ≤ z < edges[31]
    edges, d_edges = np.linspace(-cube_radius, cube_radius, num=(nbins+1), retstep=True)

    if weights is not None:
        # find which coordinates within the domain
        include = np.all(np.absolute(a) < cube_radius, axis=1)

        # find right-edge indices of non-empty bins and get which bin each coordinate belongs to
        # unique is an Mx3 array where M is the number of bins containing at least one particle and
        #   where unique[i,0], unique[i,1], and unique[i,2] give the x, y, and z indices for ith bin
        # which_bin is a 1D array giving index in unique for each coordinate in a
        unique, which_bin = np.unique(
            np.searchsorted(edges, a[include], side='right'),
            return_inverse=True,
            axis=0)

        # get weights in the domain
        w = weights[include]
        # sum the weights inside each bin
        counts = np.fromiter((w[which_bin == i].sum() for i in range(len(unique))),
                             dtype='float', count=len(unique))

    else:
        # unique is an Mx3 array where M is the number of bins containing at least one particle and
        #   where unique[i,0], unique[i,1], and unique[i,2] give the x, y, and z indices for ith bin
        # counts is a one-dimensional array of same length as unique giving count in each bin
        unique, counts = np.unique(
            np.searchsorted(edges, a[np.all(np.absolute(a) < cube_radius, axis=1)], side='right'),
            return_counts=True,
            axis=0)

    # find the coordinates of the first bin with maximal count
    coordinates_of_maximum = edges[unique[counts.argmax()]] - 0.5 * d_edges

    # we return the results
    if return_histogram:
        return (coordinates_of_maximum, edges, d_edges, unique, counts)
    else:
        return coordinates_of_maximum


# --- locate_peak_density_3D_and_plot
# same as locate_peak_density_3D, but also plot 2D histogram(s) showing the result
#   axes            an Axes to plot on or an iterable of Axes
#   squish_along    0, 1, or 2 to plot looking along x, y, or z direction; or iterable of these
#   rasterized      True or False; whether plotted histogram is raster image instead of vector
#   mark_maximum    True or False; whether to plot a marker showing location of peak density
#   n_ticks         how many tick marks to use along each axis
#
# Others as for locate_peak_density_3D. axes and squish_along should be same length if iterables.
#
def locate_peak_density_3D_and_plot(a, cube_radius, nbins, axes,
                                    squish_along=2, rasterized=True, mark_maximum=True, n_ticks=7,
                                    weights=None, return_histogram=False):
    # compute the histogram, find the peak
    hist_results = locate_peak_density_3D(a, cube_radius, nbins,
                                          weights=weights, return_histogram=True)

    # plot the 2D histogram(s)
    # is axes an Axes or an iterable?
    if hasattr(axes, '__iter__'):
        # case: iterable
        for a, sa in zip(axes, squish_along):
            _plot_results_of_locate_peak_density_3D(
                a, sa, hist_results, rasterized, mark_maximum, n_ticks)
    else:
        # case: just an Axes
        _plot_results_of_locate_peak_density_3D(
            axes, squish_along, hist_results, rasterized, mark_maximum, n_ticks)

    # return same output as locate_peak_density_3D would
    if return_histogram:
        return hist_results
    else:
        return hist_results[0]
    
# given output from locate_peak_density_3D, plot a histogram onto ax
#   ax is an Axes on which to plot
#   squish_along is 0, 1, or 2
#   hist_results is a tuple of output from locate_peak_density_3D for return_histogram=True
#   rasterized is True or False; whether plotted histogram is raster image instead of vector
#   mark_maximum is True or False; whether to indicate the peak weight density with a marker
#   n_ticks is how many tick marks to use along each axis
# doesn't return anything
cmap = copy.copy(get_cmap('plasma')) # colormap to use
cmap.set_bad(color = 'w') # changes cmap to use white background
def _plot_results_of_locate_peak_density_3D(ax, squish_along, hist_results,
                                            rasterized, mark_maximum, n_ticks):
    # extract histogram output into variables
    coordinates_of_maximum, edges, _, unique, counts = hist_results
    nbins = len(edges) - 1 # how many bins along each direction

    # get indices for the directions that are not squished
    keep_axis = tuple(i for i in range(3) if i != squish_along)

    # squish 3D histogram along the given direction to get 2D histogram
    hist2d = np.fromiter(
        (np.sum(counts[(unique[:,keep_axis[0]]==xi)
                & (unique[:,keep_axis[1]]==yi)])
         for yi in range(nbins) for xi in range(nbins)),
        dtype=counts.dtype,
        count=nbins**2
    ).reshape((nbins, nbins))
    # rows: bins along second kept axis (becomes y-axis on plot)
    # cols: bins along first kept axis (becomes x-axis on plot)

    # plot the histogram
    ax.set_xlim(edges[[0, -1]])
    ax.set_ylim(edges[[0, -1]])
    ax.set_aspect('equal')
    ax.pcolormesh(edges, edges, hist2d, norm=LogNorm(), cmap=cmap, rasterized=rasterized)
    #ax.imshow(hist2d, extent=edges[[0, -1, 0, -1]], origin='lower', norm=LogNorm(), cmap=cmap)
    # ^- imshow is an alternate raster plotter, but it requires a high figure DPI to have been set
    if mark_maximum:
        ax.scatter(coordinates_of_maximum[keep_axis[0]], coordinates_of_maximum[keep_axis[1]],
                   c='chartreuse', s=10, linewidths=0.5, marker='+')
    ticks = np.linspace(edges[0], edges[-1], n_ticks)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)

