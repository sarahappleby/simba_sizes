
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
h1_limit = 1.e6 # Msun/pc**2

rotate_galaxies = False
xlabel = 'R (kpc)'
ylabel = 'HI surface density (Msun / kpc^2)'

model = sys.argv[1]
wind = sys.argv[2]
snap = sys.argv[3]
results_dir = sys.argv[4]

results_dir += '/'+model + '_' + snap + '/'
if rotate_galaxies:
    results_dir += '/rotated_faceon'
else:
    results_dir += '/random_orientation'
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
    os.makedirs(results_dir+'/images')
    os.makedirs(results_dir+'/profiles')

data_dir = '/home/rad/data/'+model+'/'+wind+'/'
halflight_file = '/home/sapple/simba_sizes/sizes/data/halfradius_R.h5'
snapfile = data_dir+'snap_'+model+'_'+snap+'.hdf5'

sim =  caesar.load(data_dir+'Groups/'+model+'_'+snap+'.hdf5')

h = sim.simulation.hubble_constant
redshift = sim.simulation.redshift

gal_h1m = np.array([i.masses['HI'].in_units('Msun') for i in sim.galaxies])
gal_sm = np.array([i.masses['stellar'].in_units('Msun') for i in sim.galaxies])
gal_sfr = np.array([i.sfr.in_units('Msun/yr') for i in sim.galaxies])
gal_sfr[np.where(gal_sfr == 1.)[0]] = 0.
gal_ssfr = gal_sfr / gal_sm

with h5py.File(halflight_file, 'r') as f:
    gal_rad = f[model+'_'+wind+'_'+snap+'_halflight'][:] # these are in pkpc
gal_rad = np.sum(gal_rad, axis=0) / 3.

gas_pos = readsnap(snapfile, 'pos', 'gas', suppress=1, units=1) / (h*(1.+redshift)) # in kpc
gas_mass = readsnap(snapfile, 'mass', 'gas', suppress=1, units=1) / h # in Mo
gas_h1 = readsnap(snapfile, 'NeutralHydrogenAbundance', 'gas', suppress=1, units=1)

bh_pos = readsnap(snapfile, 'pos', 'bndry', suppress=1, units=1) / (h*(1.+redshift)) # in kpc

h1_radii = np.zeros(len(sim.central_galaxies))

for i in range(len(sim.central_galaxies)):
    print '\n'
    print 'Galaxy ' +str(i)

    rhalf = gal_rad[i]
    title = 'log M* = ' + str(round(gal_sm[i], 2)) + '; log(sSFR) = ' + format(round(gal_ssfr[i], 2)) + '; r200 = ' + str(round(rhalf, 2)) 
       
    glist = sim.central_galaxies[i].halo.glist
    pos = gas_pos[glist]
    mass = gas_mass[glist]
    h1 = gas_h1[glist]

    if len(sim.central_galaxies[i].bhlist) > 0.:
        pos -= bh_pos[sim.central_galaxies[i].bhlist[0]]
    else:
        print 'No black hole to center on; moving on'
        continue
        
    r  = np.sqrt(pos[:, 0]**2 +pos[:, 1]**2) 

    h1_profile, h1_radius = make_hi_profile(r, dr, h1*mass, h1_limit)
    h1_radii[i] = h1_radius[-1]

    make_image(pos[:, 0], pos[:, 1], h1*mass, results_dir+'/images/gal_'+str(i)+'.png')
    
    plt.plot(h1_radius, np.log10(h1_profile), linestyle='--', marker='.')
    plt.axhline(np.log10(h1_limit), linestyle='--', color='m')
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    plt.title(title)
    plt.xlim(0, )
    plt.savefig(results_dir+'/profiles/h1_profile_gal_'+str(i)+'.png')
    plt.clf()


with h5py.File(results_dir+model+'_'+wind+'_'+snap+'_h1_data.h5', 'a') as f:
    f.create_dataset('h1_radii', data=np.array(h1_radii))
    f.create_dataset('h1_mass', data=np.array(gal_h1m))
