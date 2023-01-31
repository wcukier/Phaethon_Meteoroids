# gererateBeta.py
# author: Wolf Cukier
# Generates a set of beta values for the geminids model.  Assumes a model run 
# total of 10,000

import pandas as pd
import numpy as np
import scipy
from scipy.interpolate import interp1d
from . import beta2mass

au = 1.495978707e11 # m/AU




s = 1.68
m1 = 1e-2


def calc_speed(m):
    """
    Assigns a speed to each particle based on its mass and a maxwellian
    distribution

    Args:
        m (ndarray, g): The masses of the particles

    Returns:
        speeds (ndarray, m/s): The speeds of the particles
    """
    v_c = 1/.4
    w=1
    v_a = v_c * np.power(m, -0.1)
    v_r = v_a.copy()

    
    for i in range(len(m)):
        done = False
        while done==False:
            x = 4*np.random.rand()
            if ((x**2*np.exp(-x**2+1))**w > np.random.rand()):

                v_a[i] = v_a[i]*x
                done = True

    return v_a

def particles(n, model, max_b=20):
    """
    Generates n particles properly spaced in beta-space with appropriate
    relative velocities for each particle

    Args:
        n (int): The number of particles
        model (int): The integer that corresponds with the model type
                            - 0 or 3: Base Model
                            - 1 or 4: Violent Creation Model
                            - 2 or 5: Cometary Creation Model
                            - 0-2: Asteroidal Composition
                            - 3-5: Young Cometary Composition
        max_b (int, optional): The maximum value for beta.  All returned
                                particles will be less than this Defaults to 20.

    Returns:
        beta (ndarray, dimentionless): The beta values of the particles
        mass (ndarray, g): The mass of the particles
        vel (ndarray, m/s): The velocities of the particles
    """
    if (max_b < .06): m0 = 1e-7;
    
    beta = np.array([])
    mass = np.array([])
    speed = np.array([])
    
    
    while beta.shape[0] < n:
        if (model % 3 == 0): b = np.linspace(0, max_b, n+1)[1:]
        else: b = np.tile(np.linspace(0, max_b, int(n/100)+1)[1:], 100)

        if (model / 3 != 1): m = beta2mass.asteroidal(b) # not the cometary comp
        else: m = beta2mass.young_comet(b) # cometary comp
        
        v = calc_speed(m)
        
        beta = np.concatenate((beta, b))
        mass = np.concatenate((mass, m))
        speed = np.concatenate((speed, v))
        
    beta = beta[:n]
    mass = mass[:n]
    speed = speed[:n]
        

    theta = np.random.rand(n) * 2 * np.pi
    phi = np.random.rand(n) * np.pi
    
    vel = np.zeros((n,3))
    vel[:,0] = speed * np.cos(theta) * np.sin(phi)
    vel[:,1] = speed * np.sin(theta) * np.sin(phi)
    vel[:,2] = speed * np.cos(phi)
    vel = vel / au * 1000
    
    return beta, mass, vel
