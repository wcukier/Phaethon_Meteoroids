# geminids/beta2mass.py
# Author: Wolf Cukier
# Converts beta values to mass in grams from asteroidal and
# young cometary models

import pandas as pd
import numpy as np
import scipy
from scipy.interpolate import interp1d
import scipy.optimize

au = 1.495978707e11 # m/AU




s = 1.68
m1 = 1e-2

# General form of beta-mass function
def fit(x,a,b):
    return a* x**-b

def asteroidal(b):
    """
    Calculates the mass the corresponds to each beta value assuming an asteroidal
    profile as defined in Wilck and Mann (1996)

    Args:
        b (ndarray): The beta values of the particles

    Returns:
        mass (ndarray, g): The corresponding masses of each of the particles
    """
    asteroid = pd.read_csv('./data/beta_asteroidal.txt', header=None) #[g], beta
    rho = 3205.64 * 1e3 #g/m^3 (calc. from Wilck and Mann)
    mass = asteroid[0]
    beta = asteroid[1]

    m = []

    params = scipy.optimize.curve_fit(fit, beta[55:], mass[55:])
    f = interp1d(beta[55:], np.log(mass[55:]), kind = "linear")

    for i in b:
        if i >= beta[55:].max(): m.append(np.nan)
        elif i >= beta.min(): m.append(np.exp(f(i)))
        else: m.append(fit(i, params[0][0], params[0][1]))

    return np.array(m)


def young_comet(b):
    """
    Calculates the mass the corresponds to each beta value assuming an young_cometary
    profile as defined in Wilck and Mann (1996)

    Args:
        b (ndarray): The beta values of the particles

    Returns:
        mass (ndarray, g): The corresponding masses of each of the particles
    """
    comet = pd.read_csv('./data/beta_young_cometary.txt', header=None) #[g], beta
    rho = 364.278 #g/m^3 (calc. from Wilck and Mann)
    mass = comet[0]
    beta = comet[1]

    m = []

    params = scipy.optimize.curve_fit(fit, beta[55:], mass[55:])
    f = interp1d(beta, np.log(mass), kind = "linear")

    for i in b:
        if i >= beta.min(): m.append(np.exp(f(i)))
        else: m.append(fit(i, params[0][0], params[0][1]))

    return np.array(m)

