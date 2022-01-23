"""
cometary_start.py
author: Wolf Cukier
Determines the initial locations of particles to be released for the cometary
model
"""
import sys
import numpy as np
from scipy.interpolate import interp1d

from .constants import *


def idx(k, ascending, decending, r_bounds, peri, n):
    """
    Returns the index of orbit the kth particle.
    """

    if (k < int(n/2)):
        return int(decending(np.flip(r_bounds)[k]))
    if ((k == int(n/2)) and (n % 2 == 1)):
        return peri
    else:
        return int(ascending(r_bounds[k-int(n/2) + (n % 2)]))

def init_loc(k, orbit, n=100):
    """
    Divides the given orbit into n points, equally spaced by radial distance
    from the sun, on either side of perihelion.  Assumes that orbit has, in this
    order, an aphelion, perihelion, and an aphelion.  Returns the kth such
    point and the relative time spent at that point.
    """

    r = np.sqrt(orbit[:,0]**2 + orbit[:,1]**2 + orbit[:,2]**2)

    peri = np.argmin(r)
    ap_1 = np.argmax(r[:peri])
    ap_2 = np.argmax(r[peri:]) + peri

    indices = np.arange(len(r))
    decending = interp1d(r[ap_1:peri+1], indices[ap_1:peri+1])
    ascending = interp1d(r[peri:ap_2+1], indices[peri:ap_2+1])

    r_bounds = np.linspace(r[peri], np.min((r[ap_1], r[ap_2])), int(n/2))

    idxk = idx(k, ascending, decending, r_bounds, peri, n)
    if (k == 0):
        t_weight = (idx(k+1, ascending, decending, r_bounds, peri, n)+idxk)/2 - ap_1
    elif (k == n-1):
        t_weight = ap_2 - (idx(k-1, ascending, decending, r_bounds, peri, n) +
                           idxk)/2
    else:
        t_weight = (idx(k+1, ascending, decending, r_bounds, peri, n) -
                    idx(k-1, ascending, decending, r_bounds, peri, n))/2

    print(t_weight, file=sys.stderr)
    return orbit[idxk], r[idxk], t_weight

def max_beta(k, orbit, n=100):
    """
    Returns the maximum beta that will survive for the kth release location
    along the given orbit assuming a total of n such locations are calculated
    using init_loc.
    """

    y0, r, t = init_loc(k, orbit, n)
    return r*AU_TO_M/(2 * PHAETHON_SEMI_MAJOR)

