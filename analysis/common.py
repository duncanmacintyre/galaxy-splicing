# common functions used in multiple other files
# DM: initial version July 2021

import copy
import warnings

import numpy as np
from matplotlib.cm import get_cmap
from matplotlib.colors import LogNorm


# --- code_time_to_Myr and code_mass_to_Msun
# for unit conversion

# given time in code units, return time in Myr
def code_time_to_Myr(t):
    return 978.028*t

# given mass in code units, return mass in M☉
def code_mass_to_Msun(m):
    return 1.e10 * m


# --- locate_peak_density_3D, locate_peak_density_3D_and_plot, and plot_2D_histogram
# for finding peak density and plotting 2D histograms

# we use this colormap in the plotting functions below
cmap = copy.copy(get_cmap('plasma'))
cmap.set_bad(color = 'w') # changes cmap to use white background

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


# same as locate_peak_density_3D, but also plot 2D histogram(s) showing the result
#   axes            an Axes to plot on or an iterable of Axes
#   squish_along    0, 1, or 2 to plot looking along x, y, or z direction; or iterable of these
#   rasterized      True or False; whether plotted histogram is raster image instead of vector
#   nticks          how many tick marks to use along each axis
#   mark_maximum    True or False; whether to plot a green plus showing location of peak density
#
# Others as for locate_peak_density_3D. axes and squish_along should be same length if iterables.
#
# Notes:
#   This funciton will be much slower than than plot_2D_histogram so the latter should be preferred
#   when mark_maximum is False and we don't need the 3D histogram that this function gives.
#
#   The maximum that this function plots is the max in the 3D space whereas the max plotted by
#   plot_2D_histogram is the max in a 2D projection. Therefore, when marking maxima, this function
#   gives a more useful/robust maximum but is much slower.
#
def locate_peak_density_3D_and_plot(a, cube_radius, nbins, axes,
                                    squish_along=2, rasterized=True, nticks=7, mark_maximum=True, 
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
                a, sa, hist_results, rasterized, mark_maximum, nticks)
    else:
        # case: just an Axes
        _plot_results_of_locate_peak_density_3D(
            axes, squish_along, hist_results, rasterized, mark_maximum, nticks)

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
#   nticks is how many tick marks to use along each axis
# doesn't return anything
def _plot_results_of_locate_peak_density_3D(ax, squish_along, hist_results,
                                            rasterized, mark_maximum, nticks):
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
    #if mark_cm:
    #    cm = np.sum(coords*weights.reshape((-1,1)), axis=1)/np.sum(weights)
    #    ax.scatter(cm[keep_axis[0]], cm[keep_axis[1]],
    #               c='red', s=6, linewidths=0.5, marker='.')
    ticks = np.linspace(edges[0], edges[-1], nticks)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)


# plot a 2D histogram "map", for example, of stellar mass
#   a                 Nx3 array a of coordinates
#   cube_radius       plot will show domain (-cube_radius, cube_radius) along each axis
#   nbins             plot will show nbins^2 boxes
#   axes              an Axes to plot on or an iterable of Axes
#   squish_along      0, 1, or 2 to plot looking along x, y, or z direction; or iterable of these
#   rasterized        True or False; whether plotted histogram is raster image instead of vector
#   nticks            how many tick marks to use along each axis
#   mark_maximum      True or False; whether to plot a green plus showing bin with most items
#   weights           None or an array of shape (N,) with weightings for points (i.e. masses)
#   return_histogram  whether to return histogram
#
# if return_histogram is False, returns nothing
#
# if return_histogram is True and mark_maximum is False, returns a length-two tuple with:
#   edges                         numpy array of bin edges, shape (nbins+1,)
#   hist                          array of histogram probability density function
# 
# If axes and squish_along are iterables, hist will be tuples of histograms.
#
# Notes:
#   This funciton will be faster than locate_peak_density_3D_and_plot and should always be
#   preferred when mark_maximum is False and we don't need the 3D histogram that the latter gives.
#
#   The maximum that this function plots is the max in the 2D projection whereas the max plotted by
#   locate_peak_density_3D_and_plot is the max in 3D space. Therefore, when marking maxima,
#   locate_peak_density_3D_and_plot gives a more useful/robust maximum but is much slower.
#
def plot_2D_histogram(a, cube_radius, nbins, axes,
                      squish_along=2, rasterized=True, nticks=7, mark_maximum=False, 
                      weights=None, return_histogram=False):
    
    # these are the bin edges
    edges = np.linspace(-cube_radius, cube_radius, nbins+1)
    d_edges = 2.*cube_radius/nbins # width of bins

    # compute and plot the 2D histogram(s)
    # is axes an Axes or an iterable?
    if hasattr(axes, '__iter__'):
        # case: iterable - we iterate over axes and squish_along
        if return_histogram:
            return (edges, tuple((_plot_2D_histogram_helper(
                                        a, edges, d_edges, nbins, ax, sa, rasterized=rasterized,
                                        nticks=nticks, mark_maximum=mark_maximum,
                                        weights=weights, return_histogram=True)
                                 for ax, sa in zip(axes, squish_along))))
        else:
            for ax, sa in zip(axes, squish_along):
                _plot_2D_histogram_helper(a, edges, d_edges, nbins, ax, sa, rasterized=rasterized,
                                          nticks=nticks, mark_maximum=mark_maximum,
                                          weights=weights, return_histogram=False)
    else:
        # case: just an Axes - we call _plot_2D_histogram_helper just once
        if return_histogram:
            return (edges, _plot_2D_histogram_helper(
                                        a, edges, d_edges, nbins, axes, squish_along=squish_along, 
                                        rasterized=rasterized, nticks=nticks,
                                        mark_maximum=mark_maximum, weights=weights, 
                                        return_histogram=True))
        else:
            _plot_2D_histogram_helper(a, edges, d_edges, nbins, axes, squish_along=squish_along, 
                                      rasterized=rasterized, nticks=nticks,
                                      mark_maximum=mark_maximum, weights=weights, 
                                      return_histogram=False)


# helper for above function - squish_along and axes can't be iterables
def _plot_2D_histogram_helper(a, edges, d_edges, nbins, ax, squish_along=2, rasterized=True, nticks=7,
                              mark_maximum=False, weights=None, return_histogram=False):
    
    # these are the coordinates along the two dimensions we keep
    x = a[:,1] if squish_along==0 else a[:,0]
    y = a[:,1] if squish_along==2 else a[:,2]
    
    # set up axes limits
    ax.set_xlim(edges[[0, -1]])
    ax.set_ylim(edges[[0, -1]])
    ax.set_aspect('equal')

    # compute and plot histogram
    if weights is not None:
        hist, _, _, _ = ax.hist2d(x, y, bins=edges, density=True, weights=weights, cmap=cmap,
                                  norm=LogNorm(), rasterized=rasterized)
    else:
        hist, _, _, _ = ax.hist2d(x, y, bins=edges, density=True, cmap=cmap,
                                  norm=LogNorm(), rasterized=rasterized)

    # find and mark the maximum
    if mark_maximum:
        index_of_maximum = np.argmax(hist) # index of maximum in flattened array
        coordinates_of_maximum = edges[list(divmod(index_of_maximum, nbins))] + d_edges
        ax.scatter(*coordinates_of_maximum, c='chartreuse', s=10, linewidths=0.5, marker='+')

    # set ticks
    ticks = np.linspace(edges[0], edges[-1], nticks)
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)

    if return_histogram:
        return hist


# ----- grab_property and empty_array
# for loading .hdf5 files (e.g. GIZMO snapshots)

default_num_metals = 11    # default number of metals
field_size = {    # how many items per field
    'Coordinates':              3,
    'Masses':                   1,
    'Velocities':               3,
    'Density':                  1,
    'InternalEnergy':           1,
    'StellarFormationTime':     1,
    'Metallicity':              default_num_metals
}

# given h5py.File f, return values for a part type and field, or an empty array if not present
#
#   num_metals   
#           how many columns in metallicity field; the number of metals tracked is num_metals-1
#   enfore_num_metals
#           if True, make sure that any returned metallicities have expected number of columns - if 
#           we load an unexpected shape, we discard it, give a RuntimeWarning, and return an array 
#           of the  expected shape where we assume zero metallicity and 25% helium
def grab_property(f, part_type, field, num_metals=None, enforce_num_metals=True):
    try:
        x = np.asarray(f['/PartType%d/%s' % (part_type, field)])
    except KeyError: # handle case where the field is not present in the hdf5 file
        if field=='Masses': # if field is 'Masses', we can try to use a mass table from the header 
            try:
                mass_per_particle = f['/Header'].attrs['MassTable'][part_type]
                n_particles = len(grab_property(f, part_type, 'Coordinates'))
                x = np.ones((n_particles,))*mass_per_particle
            except KeyError: # there was no mass table present: return empty array of correct shape
                x = empty_array('Masses')
        else: # the field is missing and isn't 'Masses': return empty array of correct shape
            x = empty_array(field)
    else: # this executes if we've finished try and a KeyError wasn't raised
        # if enforce_num_metals and we are loading metallicity, make sure desired number of metals
        if enforce_num_metals and field=='Metallicity':
            n = num_metals if num_metals is not None else default_num_metals
            if (n==1) and (len(x.shape)!=1):
                warnings.warn(
                    RuntimeWarning('The loaded metallicity had an unexpected number of metals. We ignore it and set the metallicity to zero.'))
                x = np.zeros((len(grab_property(f, part_type, 'Coordinates',
                                                num_metals=1, enforce_num_metals=False)),))
            elif (n>1) and ((len(x.shape)==1) or (x.shape[1]!=n)):
                warnings.warn(
                    RuntimeWarning('The loaded metallicity had an unexpected number of metals. We ignore it and set the metallicity to zero and He to 25%.'))
                x = np.zeros((len(grab_property(f, part_type, 'Coordinates',
                                                num_metals=n, enforce_num_metals=False)),
                              n))
                x[:,1] = 0.25

    return x

# return empty numpy array of suitable size for given field
def empty_array(field, num_metals=None):
    if (num_metals is None) or (field!='Metallicity'): # most of time
        try:
            n = field_size[field] # how many entries per particle
        except KeyError as e: # handle calls where invalid value was given for field
            valid_fields = set(field_size.keys()) # these are the valid choices for field
            if field in valid_fields: # field is valid, the KeyError is due to some other reason
                raise e 
            else: # an invalid field was given
                raise ValueError(
                    '{} is not a valid value for field. The choices are: {}'.format(
                        field, valid_fields))
    else: # field is 'Metallicity' and we need to use the custom number of metals that was given
        n = num_metals
    return np.empty((0,) if n==1 else (0, n))




