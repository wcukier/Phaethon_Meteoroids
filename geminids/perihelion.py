# perihelion.py
# Author: Wolf Cukier
# Backintegrates the orbit of Phaethon 3200 age years and finds the location
# of the perihelion of the orbit

import numpy as np
import pandas as pd
import rebound
import spiceypy as spice
from tqdm import tqdm

from .constants import *

def perihelion(age = 2000):
    spice.furnsh( 'data/meta.tm' )
    t_start = spice.str2et('1600-01-01')
    beg = spice.str2et('2018 A.D. Jan 1')
    year = spice.jyear()
    beg -= age*year
    end = 2*year

    pts_per_year = 365*24*600
    n_particles = 1
    n_years = 2


    [y0, lt] = spice.spkezr('2003200', t_start, 'ECLIPJ2000', 'NONE', 'SUN')
    y0 = spice.convrt(y0, "KM", "AU")


    n = 1



    sim = rebound.Simulation()


    sim.units = ('s', 'AU', 'Msun')
    sim.dt = -.001
    sim.add(m=1.)

    [yj, lt] = spice.spkezr("JUPITER BARYCENTER", t_start, 
                            "ECLIPJ2000", "NONE", "SUN")
    yj = spice.convrt(yj, "KM", "AU")
    sim.add(m=0.000954588, x=yj[0], y=yj[1], z=yj[2], 
            vx=yj[3], vy=yj[4], vz=yj[5])
    print("Added Jupiter...")

    [e_pos, lt] = spice.spkezr('EARTH', t_start, 'ECLIPJ2000', 'NONE', 'SUN')
    e_pos = spice.convrt(e_pos, 'KM', 'AU')
    sim.add(m=MASS_E/MASS_SUN, x=e_pos[0], y=e_pos[1], z=e_pos[2],
            vx=e_pos[3], vy=e_pos[4], vz=e_pos[5])
    print("Added Earth...")


    [mr_pos, lt] = spice.spkezr('4', t_start, 'ECLIPJ2000', 'NONE', 'SUN')
    mr_pos = spice.convrt(mr_pos, 'KM', 'AU')
    sim.add(m=MASS_MR/MASS_SUN, x=mr_pos[0], y=mr_pos[1], z=mr_pos[2], 
            vx=mr_pos[3], vy=mr_pos[4], vz=mr_pos[5])
    print("Added Mars...")



    [v_pos, lt] = spice.spkezr('VENUS', t_start, 'ECLIPJ2000', 'NONE', 'SUN')
    v_pos = spice.convrt(v_pos, 'KM', 'AU')
    sim.add(m=MASS_V/MASS_SUN, x=v_pos[0], y=v_pos[1], z=v_pos[2], 
            vx=v_pos[3], vy=v_pos[4], vz=v_pos[5])
    print("Added Venus...")



    [hg_pos, lt] = spice.spkezr('MERCURY', t_start, 'ECLIPJ2000', 'NONE', 'SUN')
    hg_pos = spice.convrt(hg_pos, 'KM', 'AU')
    sim.add(m=MASS_HG/MASS_SUN, x=hg_pos[0], y=hg_pos[1], z=hg_pos[2], 
            vx=hg_pos[3], vy=hg_pos[4], vz=hg_pos[5])
    print("Added Mercury...")


    n_active = len(sim.particles)




    sim.n_active = n_active
    sim.collision = "none"



    sim.move_to_hel()
    print("adding particles...")

    sim.add(x = y0[0], y=y0[1], z=y0[2], vx=y0[3], vy = y0[4], vz = y0[5])
    sim.move_to_com()




    print("time arrays...")

    Noutputs = n_years*pts_per_year
    year = spice.jyear()
    times = np.linspace(end-t_start, beg-t_start, int(Noutputs))


    t2 = np.linspace(0, end-t_start, int(1600 + 3000)*100)
    pos = np.zeros((n,int(Noutputs),3))

    sim.move_to_com()        
    ps = sim.particles      



    xy = np.zeros((int(Noutputs),6))
    xy2 = np.zeros((int(1600 + 3000)*100,3))
    d = np.zeros(int(Noutputs))
    print("simulating...")
    p = sim.particles[6]

    for i, time in enumerate(tqdm(t2)):
        sim.integrate(time)
        xy2[i] = [p.x, p.y, p.z]    
    for i, time in enumerate(tqdm(times)):
        sim.integrate(time)
        sim.move_to_hel()
        xy[i] = [p.x, p.y, p.z, p.vx, p.vy, p.vz]
        d[i] = spice.vnorm(xy[i][0:3])
        sim.move_to_com()

    start_date = t_start

    idxmin = np.argmin(d)
    t_start = times[idxmin] + start_date
    perihelion = xy[idxmin]

    return perihelion, t_start, xy


