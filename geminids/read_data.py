#read_data.py
#author: Wolf Cukier
#A set of functions to help with data processing
import numpy as np
from tqdm import tqdm
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.spatial import KDTree
import matplotlib.dates as mdates
import pandas as pd
import spiceypy as spice

from geminids.constants import AU_TO_M,s
from .weights import *
from .cometary_start import init_loc


RADIUS_EARTH = 384.4e6 #m
PEAK_DENSITY = 127/((59.2 + 28.7)/2) * 2.3e-2 /1e3**2/3600 # particles/m^2/sec
#              scale factor from fig   ZHR    (km->m)^2 h->s



cords = {"x":0, "y":1, "z":2}

def load(n: int, model: int, pth=None):
    """Reads in data

    Args:
        n (int): number of files to load in
        model (int): 0 for novel, 1 for vel, 2 for distr
        pth (str, optional):  The path to read the data in from.
                            Defaults to "../output/{model}"
                            where {model} is the string representation of
                            the model num

    Returns:
        points(ndarray(200000*n,5)): An array which cointains the x,y,z,weight and
                                either the beta or r value (r for distr, beta
                                otherwise)
        mass(float):                  The unnormalized mass of the stream
    """


    orbit = np.load("../data/orig_orbit.npy")
    if (not pth):  pth = "../output/" + ["novel", "vel", "distr"][model]
    points = []
    mass = 0
    for i in tqdm(range(n)):
        data = np.load(f'{pth}/particles{i}.npy')
        beta = np.load(f'{pth}/beta{i}.npy')
        size = np.load(f'{pth}/mass{i}.npy')
        size = asteroidal(beta)
        data = data.reshape(10000*10*2, 5)
        data = data[:,:3]



        if model==0:
            s = 5
            points.append(np.hstack((data, np.tile(size,2000).reshape(200000,1),
                                     np.tile(weight_novel(beta),
                                             2000).reshape(200000,1))))
            mass += np.sum(size[(~np.isnan(size * weight_novel(beta)))]
                           * weight_novel(beta[(~np.isnan(size * weight_novel(beta)))]))
        elif model==1:
            s = 5
            points.append(np.hstack((data, np.tile(size,2000).reshape(200000,1),
                            np.tile(weight_vel(beta),
                                    2000).reshape(200000,1))))
            mass += np.sum(size[(~np.isnan(size * weight_vel(beta)))]
                   * weight_vel(beta[(~np.isnan(size * weight_vel(beta)))]))
        elif model==2:
            s = 6
            __, r, t = init_loc(int(i/10), orbit)
            t = np.tile(t, 100)
            r = np.tile(r, 100)
            points.append(np.hstack((data, np.tile(size,2000).reshape(200000,1),
                    np.tile(weight_cometary(beta, r, t),2000).reshape(200000,1),
                    np.tile(r,2000).reshape(200000,1))))
            mass += np.sum(size[(~np.isnan(size * weight_cometary(beta, r, t)))]
                   * weight_cometary(beta[(~np.isnan(size * weight_vel(beta)))], r, t))

    mass *= 2000
    points = np.array(points)
    points = points.reshape(200000*n,s)
    points = points[~np.isnan(points).any(axis=1)]

    return points, mass

def load_all_data():
    """Loads all simulation data into 3 lists.  The first elements of the list
    always corresponds to the novel model, second to the vel model, third to the
    distr model.

    Returns:
        points(list of ndarray): The state vectors outputted by the three models
        orbital elements(list of ndarray): The orbital elements outputted by the three models
        masses(list of floats): The unnormalized masses for each model
    """
    n = [100, 100, 1000]

    points = []
    elements = []
    masses = []
    for i in range(3):
        try:
            points.append(np.load(f"../output/cached/points_{i}.npy"))
            masses.append(np.load(f"../output/cached/mass_{i}.npy")[0])
            print(f"Loaded ../output/cached/points_{i}.npy")
            print(f"Loaded ../output/cached/mass_{i}.npy")
        except Exception as e:
            print(e)
            print("Reloading data--This may take a few minutes")
            point, mass = load(n[i], i)
            points.append(point)
            masses.append(mass)
            np.save(f"../output/cached/points_{i}.npy", point)
            print(f"Loaded and cached file to ../output/cached/points_{i}.npy")
            np.save(f"../output/cached/mass_{i}.npy", np.array([mass]))
            print(f"Loaded and cached file to ../output/cached/mass_{i}.npy")

        try:
            elements.append(np.load(f"../output/cached/elements_{i}.npy"))
            print(f"Loaded ../output/cached/elements_{i}.npy")
        except Exception as e:
            print(e)
            print("Reloading data--This may take a few minutes")
            elements.append(load_elements(n[i], i))
            np.save(f"../output/cached/elements_{i}.npy", elements[i])
            print(f"Loaded and cached file to ../output/cached/elements_{i}.npy")

    return points, elements, masses



def get_phaethon_orbit():
    """Returns the orbit of (3200) Phaethon as an array

    Returns:
        orbit_phaethon (ndarray): Position vectors of (3200) Phaethon in ECLIP_J2000 cords
    """
    orbit = np.load("data/phaethon.npy")
    orbit = orbit.reshape((20000,5))[:,:3].copy()
    transform = spice.pxform("J2000", "ECLIPJ2000", 0)
    return np.dot(transform, orbit.T).T

def get_parker_orbit():
    """Returns the orbit of Parker Solar Probe near perihelion during orbit 4
    in ecliptic J2000 cords

    Returns:
        orbit_parker(ndarray): Position vectors of PSP
    """
    psp = np.load("data/psp.npy")
    psp = psp[:,:3].copy()
    transform = spice.pxform("J2000", "ECLIPJ2000", 0)
    return np.dot(transform, psp.T).T

def plot_column_density(points, orbit=None, psp=None, grid_size=50, plane='xy',
                        max_mass=np.inf, min_mass=-np.inf, extent=(-.5,2,-.5,2)):
    """Plots the colum density of the given points

    Args:
        points (ndarray): the points generated by the n-body simulations
        orbit (ndarray, optional): the orbit of (3200) Phaethon.
                                    Defaults to None.
        psp (ndarray, optional): the orbit of psp. Defaults to None.
        grid_size (int, optional): grid resolution. Defaults to 50.
        plane (str, optional): Which J2000 Plane to plot column densities in.
                               Must be a tring containing two of the following
                               characters in any order: 'x', 'y', 'z'. Defaults to 'xy'.
        max_mass (float, optional): The maximum limiting mass to plot.  Defaults to inf.
        min_mass (float, optional): The minimum limiting mass to plot.  Defaults to -inf.
        extent (tuple, optional): extent argument passed to plt.hexbin. Defaults to (-.5,2,-.5,2).
    """
    i = cords[plane[0]]
    j = cords[plane[1]]
    print(points.shape)
    mask = (points[:, 3] < max_mass) * (points[:, 3] > min_mass)
    print(points[mask].shape)

    x = np.linspace(extent[0], extent[1], grid_size*2)
    y = np.linspace(extent[2], extent[3], grid_size*2)
    c = np.zeros((grid_size*2)**2)
    nx, ny = np.meshgrid(x, y)

    poly = plt.hexbin(np.hstack((points[mask,i], nx.flatten())),
                      np.hstack((points[mask,j], ny.flatten())),
                      cmap = "plasma",
                      extent = extent,
                      gridsize=grid_size,
                      C = np.hstack((points[:,4], c)),
                      reduce_C_function = np.sum)
    plt.close()
    offsets = poly.get_offsets()
    C = poly.get_array()


    plt.hexbin(offsets[:,0], offsets[:,1], C=C/np.max(C),
               cmap="plasma", gridsize=grid_size)






    plt.colorbar(label = "Relative Column Density")
    if type(orbit) == np.ndarray: plt.plot(orbit[:,i], orbit[:,j],
                       c='g', label = "Phaethon 3200")
    if type(psp) == np.ndarray: plt.plot(psp[:,i], psp[:,j],
                     c='orange', label = "Parker Solar Probe")

    plt.xlabel(f"{plane[0]} (au)")
    plt.ylabel(f"{plane[1]} (au)")
    return

def view_perihelion(points, orbit, psp):
    """Plots the output near perihelion

    Args:
        points (ndarray): The simulation output
        orbit (ndarray): The orbit of (3200) phaethon
        psp (ndarray): The orbit of PSP
    Returns:
        None
    """
    peri = _find_perihelion(orbit)

    grid=50
    extent=(-.2,0.1,-.2,.1)

    x = np.linspace(extent[0], extent[1], grid*2)
    y = np.linspace(extent[2], extent[3], grid*2)
    c = np.zeros((grid*2)**2)
    nx, ny = np.meshgrid(x, y)

    poly = plt.hexbin(np.hstack((points[:,0], nx.flatten())),
                      np.hstack((points[:,1], ny.flatten())),
                      cmap = "plasma",
                      extent = extent,
                      gridsize=grid, C =
                      np.hstack((points[:,4], c)),
                      reduce_C_function = np.sum)
    plt.close()
    offsets = poly.get_offsets()
    C = poly.get_array()

    dx = (np.unique(np.sort(offsets, axis=0)[:,0])[1]
          - np.unique(np.sort(offsets, axis=0)[:,0])[0])
    dy = (np.unique(np.sort(offsets, axis=0)[:,1])[1] -
          np.unique(np.sort(offsets, axis=0)[:,1])[0])



    plt.hexbin(offsets[:,0],
               offsets[:,1],
               C=C/np.max(C),
               cmap="plasma",
               gridsize=grid,
               extent = extent)

    plt.colorbar(label = "Column Density (particles/km$^2$)")
    plt.plot(orbit[:,0],
             orbit[:,1],
             alpha = 1,
             c='g',
             label = "Phaethon 3200")

    plt.plot([orbit[peri,0]], [orbit[peri,1]], "r-o")
    plt.scatter([0], [0], marker="*", c='y', s=500)
    plt.plot(psp[:,0], psp[:,1], c='orange')

    plt.xlim(extent[0], extent[1])
    plt.ylim(extent[2], extent[3])

def _quadratic(x, a, b, c):
    return a*(x**2) + b * (x**1)  + c

def _find_perihelion(orbit):
    r = np.sqrt(orbit[:,0]**2 + orbit[:,1]**2 + orbit[:,2]**2)
    return np.argmin(r)



def fit_perihelion(points, orbit):
    """Attempts to fit a curve to the core of the stream near perihelion

    Args:
        points (ndarray): The simulation output
        orbit (ndarray): The orbit of (3200) Phaethon

    Returns:
        Curve params (Pretend its an object): a bunch of parameters that
                                            describe the fitted curve
    """
    peri = _find_perihelion(orbit)

    extent=(-.3,.2,-.3,.2)

    mask = ((points[:,0] > extent[0])
            * (points[:,0] < extent[1])
            * (points[:,1] > extent[2])
            * (points[:,1] < extent[3]))

    orbit_mask = ((orbit[:,0] > extent[0])
                  * (orbit[:,0] < extent[1])
                  * (orbit[:,1] > extent[2])
                  * (orbit[:,1] < extent[3]))

    perihelion = np.array((orbit[peri,0], orbit[peri,1]))
    p = np.sqrt(np.sum(perihelion**2))


    peri_novel = points[mask]
    peri_orbit = orbit[orbit_mask]
    r = np.sqrt(peri_novel[:,0]**2 + peri_novel[:,1]**2)
    orbit_r = np.sqrt(peri_orbit[:,0]**2 + peri_orbit[:,1]**2)


    dot = (peri_novel[:,0] * perihelion[0] + peri_novel[:,1] * perihelion[1] )

    orth = (peri_novel[:,:2]
            - np.array([dot / p**2 * perihelion[0],
                        dot / p**2 * perihelion[1]]).T)

    theta = np.arccos(dot/(p * r)) * orth[:, 0]/(np.abs(orth[:, 0])+1e-10)


    orbit_dot = (peri_orbit[:,0] * perihelion[0]
                 + peri_orbit[:,1] * perihelion[1])

    orbit_orth = (peri_orbit[:,:2]
                  - np.array([orbit_dot / p**2 * perihelion[0],
                              orbit_dot / p**2 * perihelion[1]]).T)
    orbit_theta = (np.arccos(orbit_dot/(p * orbit_r))
                   * orbit_orth[:, 0]
                   /(np.abs(orbit_orth[:, 0])+1e-10)) # 1e-10 to avoid div by 0

    f = interp1d(orbit_theta, orbit_r)

    midpoints = np.linspace(orbit_theta.min(), orbit_theta.max(), 50)
    thetas = (midpoints[:-1] + midpoints[1:])/2
    offsets = np.zeros(midpoints.shape[0] - 1)

    for i in range(1, len(midpoints)):
        mask = (theta<midpoints[i]) * (theta>midpoints[i-1])
        arr = np.array((r[mask] - f(theta[mask]),
                peri_novel[mask,4]))[:, (r[mask] - f(theta[mask])).argsort()]

        arr[1] = np.cumsum(arr[1])

        cume_dist = interp1d(arr[0], arr[1])
        rs = np.linspace(np.max((arr[0].min(), -4e6)), np.min((arr[0].max(), 2e6)), 50)
        # rs = np.linspace(-4e6, 2e6, 50)

        diff_dist = cume_dist(rs[1:])-cume_dist(rs[:-1])
        offsets[i-1] = (rs[np.argmax(diff_dist)]
                        + rs[np.argmax(diff_dist) + 1])/2


    offset_mask = (theta < orbit_theta.max()) * (theta > orbit_theta.min())
    offset_rs = r[offset_mask] - f(theta[offset_mask])

    return (r,
            theta,
            orbit_r,
            orbit_theta,
            f,
            offsets,
            thetas,
            offset_mask,
            offset_rs,
            peri_novel)

def plot_roughfit(points, peri_fit):
    r, theta, orbit_r, orbit_theta, f, offsets, thetas, __ = peri_fit
    plt.hexbin(r * np.cos(theta),
               r* np.sin(theta),
               C = points[:,4],
               reduce_C_function = np.sum,
               cmap="plasma",
               gridsize=50)

    plt.plot((offsets+f(thetas))*np.cos(thetas),
             (offsets+f(thetas))*np.sin(thetas),
             "g")

def plot_residual_fit(peri_fit, color="cyan"):
    (
        r,
        theta,
        orbit_r,
        orbit_theta,
        f,
        offsets,
        thetas,
        offset_mask,
        offset_rs,
        peri_novel) = peri_fit

    [a,b,c] , corr = curve_fit(_quadratic,
                               thetas/np.pi*180,
                               offsets*AU_TO_M/1000)

    plt.hexbin(theta[offset_mask]/np.pi*180,
               offset_rs*AU_TO_M/1000/1e5,
               cmap="plasma",
               extent=(orbit_theta.min()/np.pi*180,
                       orbit_theta.max()/np.pi*180,
                       -.03*AU_TO_M/1000/1e5,
                       .02*AU_TO_M/1000/1e5),
               gridsize=30)


    plt.scatter(thetas/np.pi*180,
                offsets*AU_TO_M/1000/1e5,
                c=color,
                s=7)

    plt.plot(thetas/np.pi*180,
             _quadratic(thetas/np.pi*180, a,b,c)/1e5,
             "k--",
             label = rf"y = {a:.2e}$x^2$ + {b:.2e}$x$ + {c:.2e}")

    plt.legend()
    plt.xlabel(r"Mean Anomonly ($^\circ$)")
    plt.ylabel(r"Radial Offset from 3200 Phaethon ($10^5$ km)")
    return corr

def plot_smoothed_fit(peri_fit):
    (r,
        theta,
        orbit_r,
        orbit_theta,
        f,
        offsets,
        thetas,
        offset_mask,
        offset_rs,
        peri_novel) = peri_fit

    [a,b,c] , corr = curve_fit(_quadratic,thetas, offsets)
    plt.hexbin(r * np.cos(theta),
               r* np.sin(theta),
               C = peri_novel[:,4],
               reduce_C_function = np.sum,
               cmap="plasma")

    plt.plot((_quadratic(thetas, a, b, c)+f(thetas))*np.cos(thetas),
             (_quadratic(thetas, a, b, c)+f(thetas))*np.sin(thetas),
             "--",
             color="darkgreen",
             linewidth=2,
             label="Fitted Stream Core")

    plt.plot(orbit_r * np.cos(orbit_theta),
             orbit_r* np.sin(orbit_theta),
             color="y",
             linewidth=3,
             label="3200 Phaethon")

    plt.xlim(0, 0.15)
    plt.ylim(-0.25, 0.25)
    plt.legend(loc=6)
    plt.xlabel("x (au)")
    plt.ylabel("y (au)")

##  Elements ###################################################################

def load_elements(n:int, model:int, pth=None):
    """
    Reads in the orbital elements of the model

    Args:
        n (int): number of files to load in
        model (int): 0 for novel, 1 for vel, 2 for distr
        pth (str, optional): Path to load data from. Defaults to "../output/{model}"
                where {model} is the string representation of the model num

    Returns:
        orbital_elements(ndarray(200000*n,8)): array that contains the following data for each
                                particle (second index):
                                    - 0: semi-major axis (au)
                                    - 1: eccentricity (dimentionless)
                                    - 2: inclination (degrees)
                                    - 3: argument of pericenter (degrees)
                                    - 4: longitude of ascending node (degrees)
                                    - 5: mass (grams)
                                    - 6: beta (dimentionless)
                                    - 7: weight (as in how many particles it reps)
    """


    orbit = np.load("../data/orig_orbit.npy")
    if (not pth):  pth = "../output/" + ["novel", "vel", "distr"][model]
    elements = np.zeros((1998,100*n,8))
    for i in tqdm(range(n)):
        data = np.load(f"{pth}/elements{i}.npy")
        beta = np.load(f"{pth}/beta{i}.npy")
        size = asteroidal(beta)
        temp_arr = np.zeros((1998,100,8))
        temp_arr[:,:,:5] = data
        temp_arr[:,:,5] = np.tile(size,1998).reshape(1998,100)
        temp_arr[:,:,6] = np.tile(beta,1998).reshape(1998,100)

        if model==0:
            temp_arr[:,:,7] =  np.tile(weight_vel(beta),1998).reshape(1998,100)
        elif model==1:
            temp_arr[:,:,7] =  np.tile(weight_vel(beta),1998).reshape(1998,100)
        elif model==2:
                __, r, t = init_loc(int(i/10), orbit)
                t = np.tile(t, 100)
                r = np.tile(r, 100)
                temp_arr[:,:,7] =  np.tile(weight_cometary(beta, r, t),1998).reshape(1998,100)

        elements[:,100*i:100*(i+1),:] = temp_arr

    elements[:, np.any(np.isnan(elements[:,:,:]), axis=(0,2)), :] = 0
    elements[:, :, 2:5] = elements[:, :, 2:5]/np.pi * 180

    return elements

def plot_elements(elements_novel, elements_vel, elements_distr,
                  save=True, pth="../figures/elements"):
    """Takes in the elements for the novel, vel, and distr and plots
        them against time

    Args:
        elements_novel (ndarray): The orbital elements from the novel model
        elements_vel (ndarray): The orbital elements from the vel model
        elements_distr (ndarray): The orbital elements from the distr model
        save (bool, optional): if True will save a file. Defaults to True.
        pth (str, optional): where to save the file. Saves to "{pth}_{i}" where
                             i is the plotted index.
                             Defaults to "../figures/elements".
    Returns:
        None
    """
    t = np.arange(1998)



    for index in range (5):
        plt.plot(t, elements_novel[:,0,index], label = "3200 Phaethon")




        mask = (np.any(elements_novel[:,:,5] > 1e-8, axis = 0)
            * np.any(elements_novel[:,:,5] > 1e-10, axis = 0)
            * elements_novel[-1,:,0] > 0)
        if index==3:
            mask2 = elements_novel[1:,mask,index] < 180
            elements_novel[1:,mask,index][mask2] += 360

        print(np.max(elements_novel[1:,mask,index]) - np.min(elements_novel[1:,mask,index]))

        plt.plot(t[1:],
                 np.sum(elements_novel[1:,mask,index]
                        * elements_novel[1:,mask,7],axis=1)
                 / np.sum(elements_novel[1:,mask,7], axis=1),
                 label = "Basic Model")

        mask = (np.any(elements_vel[:,:,5] > 1e-8, axis = 0)
                * np.any(elements_vel[:,:,5] > 1e-10, axis = 0)
                * elements_vel[-1,:,0] > 0)

        if index==3:
            mask2 = elements_vel[1:,mask,index] < 180
            elements_vel[1:,mask,index][mask2] += 360
        print(np.max(elements_vel[1:,mask,index]) - np.min(elements_vel[1:,mask,index]))

        plt.plot(t[1:],
                 np.sum(elements_vel[1:,mask,index]
                        * elements_vel[1:,mask,7], axis = 1)
                 / np.sum(elements_vel[1:,mask,7], axis=1),
                 label="Violent Creation" )

        mask = (np.any(elements_distr[:,:,5] > 1e-8, axis = 0)
                * np.any(elements_distr[:,:,5] > 1e-10, axis = 0)
                * elements_distr[-1,:,0] > 0)
        if index==3:
            mask2 = elements_distr[1:,mask,index] < 180
            elements_distr[1:,mask,index][mask2] += 360

        print(np.max(elements_distr[1:,mask,index]) - np.min(elements_distr[1:,mask,index]))

        plt.plot(t[1:],
                 np.sum(elements_distr[1:,mask,index]
                        * elements_distr[1:,mask,7], axis = 1)
                 / np.sum(elements_distr[1:,mask,7], axis=1),
                 label ="Cometary Creation" )
        if index==4: plt.legend()

        if save: plt.savefig(f"{pth}_{index}.eps")
        plt.show()
    return

### Spectrogram Plots #########################################################
def spectogram_plot(orig_points, min_y=0, max_y=1, bins=40):
    """
    Given a set of points in ECLIPJ2000, plots the points along the orbit of
    (3200) Phaethon as a function of the radial distance from Phaethon.

    Args:
        orig_points(ndarray): The set of points to be plotted
        min_y (float, optional): The minimum distance from the orbit to plot.
                                Defaults to 0.
        max_y (float, optional): The maximum distance from the orbit to plot.
                                Defaults to 1.
        bins(int, optional): The number of bins to pass to plt.hexbin()
    Returns:
        None
    """

    orbit = get_phaethon_orbit()
    diffs = np.array( np.sqrt(np.sum([(orbit[1:,i] - orbit[:-1,i]) ** 2 for i in range(3)], axis=0)) )
    arc = np.cumsum(diffs)
    r = np.sqrt(orbit[:,0]**2 + orbit[:,1]**2 + orbit[:,2]**2)
    peri = np.argmin(r)
    points = orig_points


    # Identify the normal vector to the plane and a vetctor in the plane
    # perpendicular to the perihelion vetor
    perihelion = np.array((orbit[peri,0],
                           orbit[peri,1],
                           orbit[peri, 2]))
    point_2 = np.array((orbit[peri-10000,0],
                        orbit[peri-10000,1],
                        orbit[peri-10000, 2]))
    plane = np.cross(perihelion, point_2)
    plane = plane/np.sqrt(np.sum(plane**2))
    n = np.cross(perihelion, plane)
    n = n/np.sqrt(np.sum(n**2))



    p = np.sqrt(np.sum(perihelion**2))
    r = np.sqrt(points[:,0]**2 + points[:,1]**2 + points[:,2]**2)
    orbit_r = np.sqrt(orbit[:,0]**2 + orbit[:,1]**2 + orbit[:,2]**2)



    # Assign each particle a theta value in the plane of the orbit
    ndot= np.sum([points[:,j] * n[j] for j in range(3)], axis=0)
    pdot= np.sum([points[:,j] * plane[j] for j in range(3)], axis=0)
    proj = np.array([points[:,j] - pdot * plane[j] for j in range(3)])
    r = np.sqrt(proj[0]**2 + proj[1]**2 + proj[2]**2)
    dot = (proj[0] * perihelion[0] + proj[1] * perihelion[1] + proj[2] * perihelion[2] )
    ndot[ndot >= 0] = 1
    ndot[ndot <= 0] = -1
    theta = np.arccos(dot/(p * r)) * ndot

    # """
    # Remove the following
    # """
    # grid = 35
    # x = np.linspace(-2.5, 0.5, grid*2)
    # y = np.linspace(-0.75, 0.75, grid*2)
    # c = np.zeros((grid*2)**2)
    # nx, ny = np.meshgrid(x, y)

    # poly = plt.hexbin(np.hstack((r*np.cos(theta), nx.flatten())), np.hstack((r*np.sin(theta), ny.flatten())), cmap = "plasma",
    #               extent = (-2.5,0.5,-0.75,0.75), gridsize=grid, C = np.hstack((orig_points[:,4], c)),
    #               reduce_C_function = np.sum)

    # plt.close()
    # offsets = poly.get_offsets()
    # C = poly.get_array()

    # dx = np.unique(np.sort(offsets, axis=0)[:,0])[1] - np.unique(np.sort(offsets, axis=0)[:,0])[0]
    # dy = np.unique(np.sort(offsets, axis=0)[:,1])[1] - np.unique(np.sort(offsets, axis=0)[:,1])[0]

    # V = 2*dx*dy*(au**2)/(1e3**2)

    # plt.hexbin(offsets[:,0], offsets[:,1], C=C/np.max(C), cmap="plasma", gridsize=grid)
    # plt.show()

    # """
    # End Remove
    # """


   ### Determine the theta values along the orbit of (3200) Phaethon
    ndot= np.sum([orbit[:,j] * n[j] for j in range(3)], axis=0)
    pdot= np.sum([orbit[:,j] * plane[j] for j in range(3)], axis=0)
    o_proj = np.array([orbit[:,j] - pdot * plane[j] for j in range(3)])
    orbit_r = np.sqrt(o_proj[0]**2 + o_proj[1]**2 + o_proj[2]**2)
    o_dot = (o_proj[0] * perihelion[0] + o_proj[1] * perihelion[1] + o_proj[2] * perihelion[2] )
    ndot[ndot > 0] = 1
    ndot[ndot < 0] = -1
    ndot[(ndot==0) * (orbit_r > 0.5)] = -1
    ndot[(ndot==0) * (orbit_r <= 0.5)] = 1
    arg = o_dot/(p * orbit_r)
    arg[(arg > 1) * (arg < 1.05)] = 1
    orbit_theta = np.arccos(arg) * ndot


    # Clean up the data for the interpolation maps
    start = 4796
    end = -907
    mask = ((np.abs(theta)
            > np.min(
                np.hstack((-orbit_theta[:-1][start:end][orbit_theta[:-1][start:end]<0],
                            orbit_theta[:-1][start:end][orbit_theta[:-1][start:end]>0]))))
            * (theta < np.max(orbit_theta[:-1][start:end])-.01)
            * (theta > np.min(orbit_theta[:-1][start:end])+.01))

    points = orig_points.copy()
    x_offset = interp1d(orbit_theta, orbit[:,0])
    y_offset = interp1d(orbit_theta, orbit[:,1])
    z_offset = interp1d(orbit_theta, orbit[:,2])
    r_offset = interp1d(orbit_theta, orbit_r)
    radius = interp1d(orbit_theta, orbit_r*orbit_theta/(np.abs(orbit_theta) + 1e-10))


    arc_len = interp1d(orbit_theta[:-1][start:end],  arc[start:end] - arc[peri])


    points[mask,0] -= x_offset(theta[mask])
    points[mask,1] -= y_offset(theta[mask])
    points[mask,2] -= z_offset(theta[mask])
    r2 = np.sqrt(points[mask,0]**2 + points[mask,1] ** 2 + points[mask,2]**2)


    grid = bins
    extent = (arc_len(orbit_theta[:-1][start:end]).min(),
              arc_len(orbit_theta[:-1][start:end][~np.isnan(orbit_theta[:-1][start:end])]).max(),
              min_y,
              max_y)

    poly = plt.hexbin(arc_len(theta[mask]), r2,
                    cmap = 'plasma', gridsize = grid, C = points[mask,4],
                    reduce_C_function = np.sum, extent = extent)

    plt.xlabel("Arc Length [AU]")
    plt.ylabel("Offset from 3200 Phaethon [AU]")
    return


### At earth + PSP comparasons ################################################

def generate_KDTree(points, m_cutoff, d_min, d_max):
    """Takes a set of points and returns a KDTree of those points (filtered for
    only points within d_min of d_max and of at least mass cutoff)

    Args:
        points (ndarray): The set of particles generated by the model
        m_cutoff (float, optional): The minimum mass to be included in the tree.
                                    Defaults to 1e-5.
        d_min (float):                      The minumum distance cutoff in AU
        d_max (float):                      The maximum distance cutoff in AU


    Returns:

            Points (ndarray): The filtered points
            KDTree (scipy KDTree): A KDTree of the filtered points
    """
    d = np.sqrt(points[:,0]**2 + points[:,1]**2 + points[:,2]**2)
    points = points[np.logical_not(np.logical_or(d > d_max, d < d_min))]

    points_m4 = points[:,:][points[:,3] > m_cutoff]

    return points_m4, KDTree(points_m4[~np.isnan(points_m4).any(axis=1)][:,:3])

def rate_at_earth(points, KDTrees, r=0.05, n=8000, vels=None):
    """
    Takes a list of sets of points and a list of KDTrees and returns
        a list of the rate of impacts at Earth, along with the longitude at which
        said impacts took place

    Args:
        points (list of ndarray): The filtered pointes
        KDTrees (list of KDTrees): The kd tree of points
        r (float, optional): The search radius. Defaults to 0.05.
        n (int, optional): The temporal resolution. Defaults to 8000.

    Returns:
        rates(list of ndarray): The rate of impact at each moment in time
        long(ndarray): The longitude of Earth at each point in time
        t(ndarray): The times for which these calculations were done

    """
    #        vels(ndarray): The average velocity vector of the registered particles
    #     earth_vels(ndarray): vel of earth at same time as vels


    spice.furnsh("data/meta.tm")
    t_act = np.load('data/t-3200.npy')

    pos_arr = np.zeros((n,3))
    n_particles = np.zeros((len(KDTrees), n))
    rates = np.zeros((len(KDTrees), n))
    t = t_act[np.linspace(0,19999, n, dtype=int)]
    long = np.zeros(t.shape[0])
    earth_vels = np.zeros((n, 3))
    avg_vels = np.zeros((n, 3)) * np.nan

    et = spice.str2et('2020-01-01')

    eclip2eq = spice.sxform("ECLIPJ2000", "J2000", et)[3:, 3:]

    for i,j in enumerate(tqdm(np.linspace(0,19999, n, dtype=int))):
        [pos, lt] = spice.spkezr("EARTH", t_act[j], "ECLIPJ2000", "NONE", "SUN")
        pos = spice.convrt(pos, "KM", "AU")
        pos_arr[i] = pos[:3]

        earth_loc = eclip2eq @ pos[3:]

        earth_vels[i] = earth_loc



        for k in range(len(KDTrees)):
            hits =  KDTrees[k].query_ball_point(pos[:3], r)
            n_particles[k, i] = np.sum(points[k][hits, 4])

            if len(hits) > 0 and np.any(vels != None):
                hit_pos = np.dot(vels[hits], eclip2eq)

                avg_vels[i] = np.sum(hit_pos * np.vstack([points[k][hits, 4]]*3).T, axis=0)/np.sum(points[k][hits,4])



        pos = spice.convrt(pos, "AU", "M")
        V = 4/3 * np.pi * ((r*AU_TO_M)**3)

        for k in range(len(KDTrees)):
            rates[k, i] = n_particles[k, i]/V * 35e3

        [pos_eclip, lt] = spice.spkezr("EARTH", t_act[j], "ECLIPJ2000", "NONE", "SUN")
        rad, long[i], lat = spice.reclat(pos_eclip[:3])

    return rates, long, t, avg_vels, earth_vels


def rate_at_psp(points, KDTrees, masses, r=0.05, n=8000, norm=1e15):
    """
    Takes a list of sets of points and a list of KDTrees and returns
        a list of the rate of impacts at PSP

    Args:
        points (list of ndarray): The filtered pointes
        KDTrees (list of KDTrees): The kd tree of points
        masses (ndarray): The unnormalized masses
        r (float, optional): The search radius. Defaults to 0.05.
        n (int, optional): The temporal resolution. Defaults to 8000.
        norm (float or ndarray, optional): The mass to normalize to.
                        If an array, must be the same length as points.
                        Defaults to 1e15

    Returns:
        rates(list of ndarray): The rate of impact at each moment in time
        t(ndarray): The times for which these calculations were done
    """
    p_vel = 1.1e5 # The speed of (3200) Phaethon at perihelion, m/s
    spice.furnsh("data/meta.tm")
    psp = get_parker_orbit()
    t_vel = np.load('data/t-vel.npy')
    norms = norm/np.array(masses)

    n_particles = np.zeros((len(KDTrees), n))
    rates = np.zeros((len(KDTrees), n))

    t = []
    for i in tqdm(range(int(psp.shape[0]/100))):
        t.append(t_vel[i*100])

        for k in range(len(KDTrees)):
            n_particles[k, i] = np.sum(points[k][KDTrees[k].query_ball_point(psp[i*100], r), 4])

        V = 4/3 * np.pi * ((r*AU_TO_M)**3)

        v_psp = np.sqrt(np.sum(psp[i,3:]**2))
        for k in range(len(KDTrees)):
            rates[k, i] = (n_particles[k, i]/V * norms[k]
                            * 5 * (p_vel + v_psp))


    return rates, np.array(t)


def calibrated_mass(mass, rate):
    """Returns the mass of the stream as estimated by the model.

    Args:
        mass (float, grams): The mass of the particles in the stream
        rate (ndarray, inverse seconds): The rate the particles hit earth
    Returns:
        mass (float, grams) : The estimated mass of the stream
    """
    return mass * PEAK_DENSITY/np.max(rate)

def plot_at_earth(rates, long, labels, plot_cmor=True):
    """Plots the peak-normalized fluxes of each model suplied in rates

    Args:
        rates (list of ndarray): The rates at earth for each model
        long (ndarray): The solar longitude each rate was estimated at
        labels (list of strings): The labels for each model for plotting
        plot_cmor (bool, optional): Whether or not to plot cmor data.
                                    Defaults to True.

    Returns:
        peaks(list of floats): The solar longitudes where each model peaks
        fwhm (list of floats): The full-width half-max of the rate distribution
    """
    cmor = np.genfromtxt("data/cmor.txt", delimiter = ',')
    fig = plt.figure()
    order = np.argsort(long)
    peak = np.argmax(cmor[:,1])
    d = 1
    mean = np.mean(cmor[peak-d:peak+d,1])

    style_cycle = ["solid", "dashed", "dotted", "dashdot"]
    style_index = 1

    if plot_cmor: plt.plot(cmor[:,0],
                           cmor[:,1]/mean,
                           label = 'CMOR (observed)',
                           c = 'brown')



    for i in range(len(rates)):
        plt.plot(long[order]/2/np.pi*360 + 180,
                 rates[i][order]/np.max(rates[i]),
                 label=labels[i],
                 linestyle = style_cycle[style_index])
        style_index = (style_index+1) % 4

    plt.xlim(255,267)
    plt.xlabel(r"Solar Longitude ($\lambda^\circ$)")
    plt.ylabel("Relative Flux")
    # plt.legend(loc=3)

    peaks = np.array([long[np.argmax(rate)]/2/np.pi*360+180 for rate in rates])
    left_tails = [long[:np.argmax(rate)][rate[:np.argmax(rate)]<rate[np.argmax(rate)]/2] for rate in rates]
    right_tails = [long[np.argmax(rate):][rate[np.argmax(rate):]<rate[np.argmax(rate)]/2] for rate in rates]


    fwhm = [(right_tails[i][0] - left_tails[i][-1])/2/np.pi*360 for i in range(3)]
    return peaks, fwhm


def plot_at_psp(rates, t, labels):
    """
    Plots the estimated meteoriod fluxes of each model suplied in rates at
        psp during orbit 4

    Args:
        rates (list of ndarray): The rates at earth for each model
        t (ndarray): The times the rates were calculated for
        labels (list of strings): The labels for each model for plotting


    Returns:
        None
    """
    impact = pd.read_csv('data/psp_imp_rate_orb04.txt', sep = '\s+')
    spice.furnsh("data/meta.tm")
    t3 = spice.str2et(impact["Time"])
    times = mdates.datestr2num(spice.et2utc(t3, "C", 3))
    x2 = mdates.datestr2num(spice.et2utc(t, "C", 3))

    beg  = 45
    end = 60
    for i in range(len(rates)):
            plt.plot(x2,rates[i], label=labels[i])

    plt.plot(times, impact["Rate"][:], label = 'Observed')

    plt.xlim(times[beg], times[end])
    # plt.legend(loc=4)
    plt.yscale("log")
    plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%d/%y'))
    plt.gcf().autofmt_xdate()
    plt.xlabel("Date")
    plt.ylabel(r"Impact rate (s$^{-1}$)")

def adjustment_factor(limiting_mass_used, CMOR_limit=1.8e-4, m_min=1e-9,
                      m_max=10, s=s):
    """Factor to adjust the mass by to compensate for the fact that our model
    theoretically detects a larger fraction of particles than CMOR would

    Args:
        limiting_mass_used (float, grams): The limiting mass used to simulate
                                            collisions with earth
        CMOR_limit (float, grams, optional): The minimal mass CMOR can observe.
                                            Defaults to 1.8e-4.
        m_min (float, grams, optional): The minimum mass of particle used in the
                                        simulation. Defaults to 1e-9
        m_max (float, grams, optional): The maximum mass of particle used in the
                                        simulation. Defaults to 10
        s (float, dimentionless, optional): The power law index for the geminids
                                            Defaults to 1.68

    Returns:
        adjustment_factor (float, dimentionless)
    """
    # frac_registed = 1 - ((limiting_mass_used**(1-s) - m_min**(1-s))
    #                  / (m_max**(1-s) - m_min**(1-s)))

    # frac_observed = 1 - ((CMOR_limit**(1-s) - m_min**(1-s))
    #                 / (m_max**(1-s) - m_min**(1-s)))

    return (limiting_mass_used**(1-s)- m_max**(1-s))/(CMOR_limit**(1-s))