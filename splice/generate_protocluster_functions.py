# These functions are used in generate_protocluster.py. The reason for keeping 
# these in a seperate file is so that they can easily be imported for testing.

from os import listdir
import numpy as np

##### CONSTANTS #####

r_maximum = 109.28   # we use uniform distribution for distance from centre - this is the maximum, in kpc
v_std = 571.83       # standard deviation for galaxy speeds, in km/s - defines normal distribution
v_maximum = 1681.9   # maximum allowed galaxy speed, in km/s


##### FUNCTIONS #####

def rand_rotation_matrix(rng=None, deflection=1.0, randnums=None):
    """
    Creates a random rotation matrix.
    
    deflection: the magnitude of the rotation. For 0, no rotation; for 1, 
    competely random rotation. Small deflection => small perturbation.
    randnums: 3 random numbers in the range [0, 1]. If `None`, 
    they will be auto-generated.
    """
    # from http://www.realtimerendering.com/resources/GraphicsGems/
    #gemsiii/rand_rotation.c
    
    if randnums is None:
        randnums = rng.uniform(size = (3,))
        
    theta, phi, z = randnums
    
    theta = theta * 2.0 * deflection * np.pi  # Rotation about the pole (Z).
    phi = phi * 2.0 * np.pi  # For direction of pole deflection.
    z = z * 2.0 * deflection  # For magnitude of pole deflection.
    
    # Compute a vector V used for distributing points over the sphere
    # via the reflection I - V Transpose(V).  This formulation of V
    # will guarantee that if x[1] and x[2] are uniformly distributed,
    # the reflected points will be uniform on the sphere.  Note that V
    # has length sqrt(2) to eliminate the 2 in the Householder matrix.
    
    r = np.sqrt(z)
    V = (
            np.sin(phi) * r,
            np.cos(phi) * r,
            np.sqrt(2.0 - z)
        )
    
    st = np.sin(theta)
    ct = np.cos(theta)
    
    R = np.array(((ct, st, 0), (-st, ct, 0), (0, 0, 1)))
    
    # Construct the rotation matrix  ( V Transpose(V) - I ) R.
    
    M = (np.outer(V, V) - np.eye(3)).dot(R)
    return M


# return a length-3 unit vector in a random direction (all directions equally likely)
#       rng is an instance of numpy.random.Generator
def rand_unit_vector(rng):
    phi = rng.uniform(0, 2.*np.pi)
    costheta = rng.uniform(-1, 1)
    sintheta = np.sin(np.arccos(costheta))
    return np.array((sintheta * np.cos(phi), sintheta * np.sin(phi), costheta))


# GENERATOR: sample half-normal distribution with mode at 0 and yield only values <= maximum
#       rng is an instance of numpy.random.Generator
def sample_half_gaussian_within_bounds(rng, std, maximum):
    while True:
        mag = np.abs(rng.normal(0, std, 1))
        if mag <= maximum:
            yield mag


# GENERATOR: yield samples from a uniform distribution between a and b
#       rng is an instance of numpy.random.Generator
def sample_uniform(rng, a=0., b=1.):
    while True:
        yield rng.uniform(low=a, high=b)


# given Nx3 array, apply rotation matrix to each row then translate by vector offset
def set_new_coords(coords_vec, rotation, translation):
    for i, vec in enumerate(coords_vec):
        coords_vec[i, :] = np.dot(rotation, vec)
    return coords_vec + np.asarray(translation).reshape((1, 3))


# return how many files inside a given directory have names ending with suffix
def count_files_inside_with_suffix(path, suffix):
    return sum(1 for _ in filter(lambda s: s.endswith(suffix), listdir(path)))



