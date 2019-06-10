
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

rotate_galaxies = False
mass_bins = [10., 10.5, 11.]
bin_labels = ['10.0 - 10.5', '10.5 - 11.0', '> 11.0']
xlabel = 'R (kpc)'
ylabel = 'HI frac'

model = sys.argv[1]
wind = sys.argv[2]
snap = sys.argv[3]
sample_file = sys.argv[4]
results_dir = sys.argv[5]

if 'gv' in sample_file:
    selection = 'green_valley'
elif 'sf' in sample_file:
    selection = 'star_forming'

results_dir += '/'+model + '_' + snap + '/' + selection
if rotate_galaxies:
    results_dir += '/rotated_faceon'
else:
    results_dir += '/random_orientation'
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
    os.makedirs(results_dir+'/images')
    os.makedirs(results_dir+'/profiles')

with h5py.File(sample_file, 'r') as f:
    try:
        gal_ids = f[model+'_'+snap].value
    except KeyError:
        print 'Need to identify galaxies; run gv_sample.py first'

data_dir = '/home/rad/data/'+model+'/'+wind+'/'
halflight_file = '/home/sapple/simba_sizes/sizes/data/halfradius_R.h5'
snapfile = data_dir+'snap_'+model+'_'+snap+'.hdf5'

sim =  caesar.load(data_dir+'Groups/'+model+'_'+snap+'.hdf5')

h = sim.simulation.hubble_constant
redshift = sim.simulation.redshift

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

for m in range(len(mass_bins)):
    print '\n'
    print 'Looking at mass bin ' + bin_labels[m]
    if m != 2:
        sm_mask = (gal_sm > mass_bins[m]) & (gal_sm < mass_bins[m+1])
    else:
        sm_mask = gal_sm > mass_bins[m]

    gal_sm_use = gal_sm[sm_mask]
    gal_rad_use = gal_rad[sm_mask]
    gal_ssfr_use = gal_ssfr[sm_mask]
    gal_ids_use = gal_ids[sm_mask]

    no_gals[m] = len(gal_ids_use)
    print str(no_gals[m]) + ' galaxies in bin'
    print '\n'

    hi_radii = np.zeros(len(gal_ids_use))
    hi_profiles = []

    for i in range(len(gal_ids_use)):
        if not sim.galaxies[i].central:
            continue
        else:
            print '\n'
            print 'Galaxy ' +str(gal_ids_use[i])

            r200 = sim.galaxies[i].radii['r200'].in_units('kpc')
            title = 'log M* = ' + str(round(gal_sm_use[i], 2)) + '; log(sSFR) = ' + format(round(gal_ssfr_use[i], 2)) + '; r200 = ' + str(round(r200, 2)) 
        
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

            hi_profile, hi_radius = make_hi_profile(r, dr, hi*mass, hi_limit)
            hi_radii = hi_radius[-1]
            hi_profiles.append(hi_profile)
            plot_profile(hi_radius, hi_profile, results_dir+'profiles/hi_profile_gal_'+str(gal_ids_use[i])+'.png', xlabel=xlabel, ylabel=ylabel, title=title)

    hi_profiles = np.array(hi_profiles)

    with h5py.File(results_dir+'mask_'+str(m)+'_hi_profile.h5', 'a') as f:
        f.create_dataset('h1_profile', data=np.array(hi_profiles))
        f.create_dataset('h1_radii', data=np.array(hi_radii))
