# common functions used in multiple other files
# DM: initial version July 2021

import numpy as np

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

