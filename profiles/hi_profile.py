import h5py
import sys
import numpy as np
import os
import gc
import matplotlib.pyplot as plt

sys.path.append('/home/sapple/tools/')
from projection import *

from readgadget import readsnap
import caesar
from astropy.cosmology import FlatLambdaCDM

from profile_methods import *

dr = 1. # kpc
hi_limit = 1. # Msun/kpc**2

model = sys.argv[1]
wind = sys.argv[2]
snap = sys.argv[3]
results_dir = sys.argv[5]

data_dir = '/home/rad/data/'+model+'/'+wind+'/'
halflight_file = '/home/sapple/simba_sizes/sizes/data/halfradius.h5'
snapfile = data_dir+'snap_'+model+'_'+snap+'.hdf5'

sim =  caesar.load(data_dir+'Groups/'+model+'_'+snap+'.hdf5')

h = sim.simulation.hubble_constant
redshift = sim.simulation.redshift
cosmo = FlatLambdaCDM(H0=100*h, Om0=sim.simulation.omega_matter, Ob0=sim.simulation.omega_baryon, Tcmb0=2.73)
thubble = cosmo.age(redshift).value # in Gyr

gal_cent = np.array([i.central for i in sim.galaxies])
gal_sm = np.array([i.masses['stellar'].in_units('Msun') for i in sim.galaxies])
gal_sfr = np.array([i.sfr.in_units('Msun/yr') for i in sim.galaxies])
gal_sfr[np.where(gal_sfr == 1.)[0]] = 0.
gal_ssfr = gal_sfr / gal_sm

with h5py.File(halflight_file, 'r') as f:
    gal_rad = f[model+'_'+wind+'_'+snap+'_halflight'][:] # these are in pkpc
gal_rad = np.sum(gal_rad, axis=0) / 3.

gas_pos = readsnap(snapfile, 'pos', 'gas', suppress=1, units=1) / (h*(1.+redshift)) # in kpc
gas_mass = readsnap(snapfile, 'mass', 'gas', suppress=1, units=1) / h # in Mo
gas_vels = readsnap(snapfile, 'vel', 'gas', suppress=1, units=0) # in km/s
gas_h1 = readsnap(snapfile, 'NeutralHydrogenAbundance', 'gas', suppress=1, units=1)

bh_pos = readsnap(snapfile, 'pos', 'bndry', suppress=1, units=1) / (h*(1.+redshift)) # in kpc

for i in range(len(sim.galaxies)):
    if not sim.galaxies[i].central:
        continue
    else:
        r200 = sim.galaxies[i].radii['r200'].in_units('kpc')
        
        rhalf = gal_rad[i]
       
        glist = sim.galaxies[i].halo.glist
        pos = gas_pos[glist]
        vel = gas_vels[glist]
        mass = gas_mass[glist]
        h1 = gas_h1[glist]
        
        if len(sim.galaxies[i].bhlist) > 0.:
            pos -= bh_pos[sim.galaxies[i].bhlist[0]]
        else:
            print 'No black hole to center on; moving on'
            continue
        
        r  = np.sqrt(pos[:, 0]**2 +pos[:, 1]**2)

        profile = []
        radius = []
        j = 0
        stop = False
        while not stop:
            mask = (r > j*dr) * (r < (j+1)*dr)
            profile.append(np.sum(hi*mass[mask]))
            radius.append(j*dr + 0.5*dr)
            if j == 0:
                profile[-1] /= np.pi*(dr**2)
            else:
                profile[-1] /= np.pi* ( (dr*(j+1))**2 - (dr*j)**2)
            if profile[-1] <= hi_limit:
                stop = True
            else:
                j += 1

