import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.integrate import odeint
from scipy.interpolate import interp1d
from tqdm import tqdm
import scipy
import spiceypy as spice
import rebound
import reboundx
import os
import sys

GRAVITATIONAL_CONSTANT = 6.67430e-11 #m^3/kg/s^2
MASS_SUN = 1.98847e30 #kg
AU_TO_M = 1.495978707e11 # m/AU
DAY_TO_SEC = 86400 # sec/day
L_SUN = 3.828e26 # W
SPEED_OF_LIGHT = 299792458 #m/s
MASS_J = 1.89813e27 #kg
MASS_E = 5.97237e24 #kg
MASS_MR = 6.4171e23 #kg
MASS_V = 4.8675e24 #kg
MASS_HG = 3.3011e23 #kg
year = spice.jyear() # Year [s]
au = AU_TO_M # Astronomical Unit [m]

from generate import particles
from perihelion import perihelion

if (__name__ == "__main__"):
        k = int(sys.argv[1])
        model = int(sys.argv[2])
        dr = os.getcwd()
        if   (model == 0): subdir = "novel"
        elif (model == 1): subdir = "vel"
        elif (model == 2): subdir = "distr"
        else:
             print("second argument must be 0,1, or 2")
             raise
     
        print(k)
        n_years = 2000
        window = 2
        pts_per_year = 10000
        n_particles = 10

        y0, t_start, orbit = perihelion()


        for l in range(10):
                print(f"{k}:{l}")
                beta, mass, vel = particles(n_particles, "asteroidal", max_b = 0.052)


                n = n_particles
                b = np.concatenate(([0,0],beta))



                sim = rebound.Simulation()
                sim.tstep = .0001

                
                sim.units = ('s', 'AU', 'Msun')
                sim.add(m=1.)

                [yj, lt] = spice.spkezr("JUPITER BARYCENTER", t_start, "J2000", "NONE", "SUN")
                yj = spice.convrt(yj, "KM", "AU")
                sim.add(m=0.000954588, x=yj[0], y=yj[1], z=yj[2], vx=yj[3], vy=yj[4], vz=yj[5])

                [e_pos, lt] = spice.spkezr('EARTH', t_start, 'J2000', 'NONE', 'SUN')
                e_pos = spice.convrt(e_pos, 'KM', 'AU')
                sim.add(m=MASS_E/MASS_SUN, x=e_pos[0], y=e_pos[1], z=e_pos[2], vx=e_pos[3], vy=e_pos[4], vz=e_pos[5])

        
                [mr_pos, lt] = spice.spkezr('4', t_start, 'J2000', 'NONE', 'SUN')
                mr_pos = spice.convrt(mr_pos, 'KM', 'AU')
                sim.add(m=MASS_MR/MASS_SUN, x=mr_pos[0], y=mr_pos[1], z=mr_pos[2], vx=mr_pos[3], vy=mr_pos[4], vz=mr_pos[5])



                [v_pos, lt] = spice.spkezr('VENUS', t_start, 'J2000', 'NONE', 'SUN')
                v_pos = spice.convrt(v_pos, 'KM', 'AU')
                sim.add(m=MASS_V/MASS_SUN, x=v_pos[0], y=v_pos[1], z=v_pos[2], vx=v_pos[3], vy=v_pos[4], vz=v_pos[5])



                [hg_pos, lt] = spice.spkezr('MERCURY', t_start, 'J2000', 'NONE', 'SUN')
                hg_pos = spice.convrt(hg_pos, 'KM', 'AU')
                sim.add(m=MASS_HG/MASS_SUN, x=hg_pos[0], y=hg_pos[1], z=hg_pos[2], vx=hg_pos[3], vy=hg_pos[4], vz=hg_pos[5])


                n_active = len(sim.particles)




                sim.n_active = n_active
                sim.Nactive = n_active
                sim.collision = "none"



                sim.move_to_hel()
                for i in range(n):
                        if (model == 2): y0 = orbit[10407*(i+10*l)+105*k]
                        if (model != 1): vel[i] = [0,0,0]
                        sim.add(x = y0[0]+1e-10*np.random.rand(), y=y0[1], z=y0[2], vx=y0[3]+vel[i,0], vy = y0[4]+vel[i,1], vz = y0[5]+vel[i,2])
                sim.move_to_com()

                rebx = reboundx.Extras(sim)
                rf = rebx.load_force("radiation_forces")
                rebx.add_force(rf)
                rf.params["c"] = 0.002004



                Noutputs = window*pts_per_year
                year = spice.jyear()
                times = np.linspace(int((n_years-window)*year), int(n_years*year), int(Noutputs))

                sim.move_to_com()        
                ps = sim.particles      

                for i in range(n):
                        ps[i+n_active].params["beta"] = beta[i]

                xy = np.zeros((int(Noutputs), n, 5))

                np.save(f"{dr}/output/{subdir}/beta{k}.npy", beta)
                np.save(f"{dr}/output/{subdir}/mass{k}.npy", mass)

                for i in tqdm(range(int((n_years-window)))):
                        for p in sim.particles[n_active:]:
                            h = p.hash
                            if sim.particles[0] ** p < .01:
                                sim.remove(hash=h)
                                n_escaped +=1
                                print(f"Number escaped: {n_escaped}")
                
                escaped = True   
                while(escaped == True):
                        try:
                            sim.integrate(int(i*year))

                            escaped = False
                        except:
                            n_escaped += 1
                            print(f"Number escaped: {n_escaped}")
                            for j in range(sim.N):
                                p = sim.particles[j]
                                d2 = p.x*p.x + p.y*p.y + p.z*p.z
                                if d2>sim.exit_max_distance**2:
                                    index=p.hash # cache index rather than remove here since our loop would go beyond end of particles array
                            print(index)
                            try: sim.remove(hash=index)
                            except: escaped = False 
                        
                sim.save(f"{dr}/output/{subdir}/sim{k}.bin")
                for i, time in enumerate(tqdm(times)):

                    escaped = True
                    while(escaped == True):
                        try:
                            sim.integrate(time)  
                            escaped = False
                        except:
                            n_escaped += 1
                            print(f"Number escaped: {n_escaped}")
                        for j in range(sim.N):
                            p = sim.particles[j]
                            d2 = p.x*p.x + p.y*p.y + p.z*p.z
                            if d2>sim.exit_max_distance**2:
                                index=p.hash # cache index rather than remove here since our loop would go beyond end of particles array
                        print(index)
                        try: sim.remove(hash=index)
                        except: escaped = False
                            

                    for j in range(n):
                        try:
                            p = sim.particles[f"{j}"]
                            o = p.calculate_orbit(primary = ps[0])
                            xy[i][j] = [p.x, p.y, p.z, o.a, o.e]
                        except:
                            xy[i][j] = [np.nan, np.nan, np.nan, np.nan, np.nan]
            #                 outfile.write(f"{p.x}, {p.y}, {p.z}, {o.a}, {o.e} \n")
                np.save(f"{dr}/output/{subdir}/particles{k}.npy", xy)
