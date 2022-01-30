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
def cume_mass(m, m0=m0):
    return (m**(1-s) - m0**(1-s))/(m1**(1-s) - m0**(1-s))


def weight_novel(b):

    arr =  cume_mass(asteroidal(b)) - cume_mass(asteroidal(b + .052/10000))
    arr[asteroidal(b) > 10] = 0
    arr[asteroidal(b) < 1e-9] = 0

    return arr / asteroidal(b)


def weight_vel(b):
    arr =  cume_mass(asteroidal(b)) - cume_mass(asteroidal(b + .052/100))
    arr[asteroidal(b) > 10] = 0
    arr[asteroidal(b) < 1e-9] = 0

    return arr/100 /asteroidal(b)

def weight_cometary(b, r, t):

    mask = b < .8
    arr = b.copy()
    arr[mask] =  cume_mass(asteroidal(b[mask]), m0 = 1e-16) - cume_mass(
        asteroidal(b[mask] + max_beta_r(r[mask])/100), m0=1e-16)
    arr[mask][asteroidal(b[mask]) > 10] = 0
    arr[mask][asteroidal(b[mask]) < 1e-16] = 0
    arr[~mask] = 0
    return np.abs(arr/100) / (r ** 2) * t

def weight_novelc(b):

    arr =  cume_mass(young_comet(b)) - cume_mass(young_comet(b + .052/10000))
    arr[young_comet(b) > 10] = 0
    arr[young_comet(b) < 1e-9] = 0

    return arr / young_comet(b)