# weights.py
# author: Wolf Cukier
# A set of functions that take the beta of a particle and returns the
# weighting for the particle

from .beta2mass import asteroidal as asteroidal
from .beta2mass import young_comet as young_comet
import numpy as np
from .cometary_start import max_beta_r as max_beta_r

s = 1.68
m0 = 1e-9
m1 = 10

def cume_mass(m, m0=m0, m1=m1, s=s):
    """
    Returns the value of the cumulative mass distribution function for a power
    law distribution with power law index s at m

    Args:
        m (float or array-like): The mass(es) of the particles
        m0 (float, optional): The mass that should map to 0 on the distribution
                                function. Defaults to m0.
        m1 (float, optional): The mass that should map to 1 on the distribution
                                function. Defaults to m1.
        s (float, optional): The powerlaw index of the distribution function.
                            Defaults to 1.68.

    Returns:
        y (float): The value of the cumulative distribution function evaluated
                    at m
    """
    return (m**(1-s) - m0**(1-s))/(m1**(1-s) - m0**(1-s))


def weight_novel(b):
    """
    Returns the relative weight of the particles for the base model

    Args:
        b (ndarray): The beta values of the particles

    Returns:
        weights (ndarray): The relative weights of the particles
    """

    arr =  cume_mass(asteroidal(b)) - cume_mass(asteroidal(b + .052/10000))
    arr[asteroidal(b) > 10] = 0
    arr[asteroidal(b) < 1e-9] = 0

    return arr


def weight_vel(b):
    """
    Returns the relative weight of the particles for the violent creation model

    Args:
        b (ndarray): The beta values of the particles

    Returns:
        weights (ndarray): The relative weights of the particles
    """
    arr =  cume_mass(asteroidal(b)) - cume_mass(asteroidal(b + .052/100))
    arr[asteroidal(b) > 10] = 0
    arr[asteroidal(b) < 1e-9] = 0

    return arr/100

def weight_cometary(b, r, t):
    """
    Returns the relative weight of the particles for the cometery creation model

    Args:
        b (ndarray): The beta values of the particles
        r (ndarray): The radial distances from the sun the particles were released at
        t (ndarray): The amount of time the particle spent being represented by this release point

    Returns:
        weights (ndarray): The relative weights of the particles
    """
    mask = b < .8
    arr = b.copy()
    arr[mask] =  cume_mass(asteroidal(b[mask]), m0 = 1e-16) - cume_mass(
        asteroidal(b[mask] + max_beta_r(r[mask])/100), m0=1e-16)
    arr[mask][asteroidal(b[mask]) > 10] = 0
    arr[mask][asteroidal(b[mask]) < 1e-16] = 0
    arr[~mask] = 0
    return np.abs(arr/100) / (r ** 4) * t

def weight_novelc(b):
    """
    Returns the relative weight of the particles for the base model but assuming
    a young cometary composition

    Args:
        b (ndarray): The beta values of the particles

    Returns:
        weights (ndarray): The relative weights of the particles
    """
    arr =  cume_mass(young_comet(b)) - cume_mass(young_comet(b + .052/10000))
    arr[young_comet(b) > 10] = 0
    arr[young_comet(b) < 1e-9] = 0

    return arr