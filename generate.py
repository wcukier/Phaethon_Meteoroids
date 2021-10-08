# gererate.py

import pandas as pd
import numpy as np
import scipy
from scipy.interpolate import interp1d

au = 1.495978707e11 # m/AU




s = 1.68
m1 = 1e-2

def inv_mass_distr(x, m0 = 1e-9):
    c = (1-s)/(m1**(1-s) - m0 ** (1-s) )
    return ((1-s)/c*x + m0**(1-s))**(1/(1-s))

def gen_mass(n, m0 = 1e-9):
    x = np.random.rand(n)
    return inv_mass_distr(x, m0)


# General form of beta-mass function
def fit(x,a,b):
    return a* x**b

def asteroidal(m):
    asteroid = pd.read_csv('data/beta_asteroidal.txt', header=None) #[g], beta
    rho = 3205.64 * 1e3 #g/m^3 (calc. from Wilck and Mann) 
    mass = asteroid[0]
    beta = asteroid[1]

    params = scipy.optimize.curve_fit(fit, mass[55:], beta[55:])
    f = interp1d(np.log(mass),beta, kind = "linear")

    m = np.array([i for i in m if mass.min() < i < 4/3*np.pi*rho*(1e-2**2)])
    s = (m / (4/3 * np.pi * rho)) ** (1/3)
    b = []

    for i in m:
        if i <= mass.max(): b.append(f(np.log(i)))
        else: b.append(fit(i, params[0][0], params[0][1]))

    v_c = 1/.4
    w=1
    v_a = v_c * np.power(m, 0.1)
    v_r = v_a.copy()

    
    for i in range(len(m)):
        done = False
        while done==False:
            x = 4*np.random.rand()
            if ((x**2*np.exp(-x**2+1))**w > np.random.rand()):

                v_a[i] = v_a[i]*x
                done = True
    


    return np.array(b), m, v_r


def particles(n, kind, max_b=20):
    if (max_b < .06): m0 = 1e-9;
    else: m0 = 1e-16;
    
    beta = np.array([])
    mass = np.array([])
    speed = np.array([])
    
    
    while beta.shape[0] < n:
        m = gen_mass(n - beta.shape[0]) #grams

        b, m, v = asteroidal(m)

        m = m[b<= max_b]
        v = v[b<= max_b]
        b = b[b<= max_b]


        beta = np.concatenate((beta, b))
        mass = np.concatenate((mass, m))
        speed = np.concatenate((speed, v))
        


    theta = np.random.rand(n) * 2 * np.pi
    phi = np.random.rand(n) * np.pi
    
    vel = np.zeros((n,3))
    vel[:,0] = speed * np.cos(theta) * np.sin(phi)
    vel[:,1] = speed * np.sin(theta) * np.sin(phi)
    vel[:,2] = speed * np.cos(phi)
    vel = vel / au * 1000
    
    return beta, mass, vel
