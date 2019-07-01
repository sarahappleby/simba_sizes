# first choose some galaxies

# make mass cuts, 10.0 - 10.5, 10.5 - 11.0, > 11.
# at each mass take a star forming galaxy and a quenched galaxy if these exist

# for each of these galaxies:
# make image of h1, h2, stars
# make profiles of h1 mass, h2 mass, ssfr, temp, star age

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

sys.path.append('..')
from profile_methods import *

dr = 1. # kpc
factor = 5.
ssfr_limit = -10.

mass_bins = [10., 10.5, 11.]
bin_labels = ['10.0 - 10.5', '10.5 - 11.0', '> 11.0']

model = sys.argv[1]
wind = sys.argv[2]
snap = sys.argv[3]
results_dir = sys.argv[4]

data_dir = '/home/rad/data/'+model+'/'+wind+'/'
halflight_file = '/home/sapple/simba_sizes/sizes/data/halfradius_R.h5'
snapfile = data_dir+'snap_'+model+'_'+snap+'.hdf5'

sim =  caesar.load(data_dir+'Groups/'+model+'_'+snap+'.hdf5')

h = sim.simulation.hubble_constant
redshift = sim.simulation.redshift

gal_ids = np.arange(len(sim.central_galaxies))
gal_rad = np.array([i.radii['stellar_half_mass'].in_units('Msun') for i in sim.central_galaxies])
gal_sm = np.array([i.masses['stellar'].in_units('Msun') for i in sim.central_galaxies])
gal_sfr = np.array([i.sfr.in_units('Msun/yr') for i in sim.central_galaxies])
gal_ssfr = gal_sfr / gal_sm

gal_sm = np.log10(gal_sm)
gal_ssfr = np.log10(gal_ssfr)

sf_mask = gal_ssfr > ssfr_limit

gas_pos = readsnap(snapfile, 'pos', 'gas', suppress=1, units=1) / (h*(1.+redshift)) # in kpc
gas_mass = readsnap(snapfile, 'mass', 'gas', suppress=1, units=1) / h # in Mo
gas_sfr = readsnap(snapfile, 'sfr', 'gas', suppress=1, units=1) # in Mo/yr
gas_h1 = readsnap(snapfile, 'NeutralHydrogenAbundance', 'gas', suppress=1, units=1)
gas_h2 = readsnap(snapfile, 'fh2', 'gas', suppress=1, units=1)
gas_temp = readsnap(snapfile, 'u', 'gas', suppress=1, units=1)

bh_pos = readsnap(snapfile, 'pos', 'bndry', suppress=1, units=1) / (h*(1.+redshift)) # in kpc

gals = []

for j, b in enumerate(bin_labels):
    print '\n'
    print 'Looking at mass bin ' + b

    if j != len(mass_bins) - 1:
        sm_mask = (gal_sm > mass_bins[j]) & (gal_sm < mass_bins[j+1])
    else:
        sm_mask = gal_sm > mass_bins[j]

    # star forming:

    gal_ids_use = gal_ids[sm_mask*sf_mask]
    choose = np.random.randint(len(gal_ids_use))
    gals.append(gal_ids_use[choose])

    # quenched galaxies:

    gal_ids_use = gal_ids[sm_mask*np.invert(sf_mask)]
    choose = np.random.randint(len(gal_ids_use))
    gals.append(gal_ids_use[choose])

for i in gals:
    print 'Galaxy ' +str(i)
    slist = sim.central_galaxies[i].halo.slist
    glist = sim.central_galaxies[i].halo.glist
    rhalf = gal_rad_use[i]
    title = 'log M* = ' + str(round(gal_sm[i], 2)) + '; log(sSFR) = ' + format(round(gal_ssfr[i], 2))
    print title
    
    n = int(round(rhalf*factor / dr))
    rplot = np.arange(0., dr*n, dr) + (dr*0.5)

    pos = star_pos[slist]
    mass = star_mass[slist]
    tform = star_tform[slist]
    ages = tage(cosmo, thubble, tform) * 1000. # get star ages in Myr

    if len(sim.central_galaxies[i].bhlist) > 0.:
        pos -= bh_pos[sim.central_galaxies[i].bhlist[0]]
    else:
        pos -= center_of_quantity(pos, mass)
        print 'No black holes to center on, centering on stars'

    r = np.sqrt(pos[:, 0]**2 +pos[:, 1]**2)

    sm_prof = real_profile(n, dr, r, mass)
    age_prof = real_profile(n, dr, r, ages)
    plot_profile(rplot, sm_prof, results_dir+'/sm_profile_gal_'+str(i)+'.png', 'Mass surface density', title=title)
    plot_profile(rplot, age_prof, results_dir+'/age_profile_gal_'+str(i)+'.png', 'Age surface density', title=title)

    pos = gas_pos[glist]
    mass = gas_mass[glist]
    sfr = gas_sfr[glist]
    h1 = gas_h1[glist]
    h2 = gas_h2[glist]
    temp = gas_temp[glist]

    if len(sim.galaxies[gal_ids_use[i]].bhlist) > 0.:
        pos -= bh_pos[sim.galaxies[gal_ids_use[i]].bhlist[0]]
    else:
        pos -= center_of_quantity(pos, mass)
        print 'No black holes to center on, centering on stars'

    r = np.sqrt(pos[:, 0]**2 +pos[:, 1]**2)

    gas_m_prof = real_profile(n, dr, r, mass)
    sfr_prof = real_profile(n, dr, r, sfr)
    h1_prof = real_profile(n, dr, r, h1*mass)
    h2_prof = real_profile(n, dr, r, h2*mass)
    temp_prof = real_profile(n, dr, r, temp*mass)
    temp_prof / gas_m_prof
    ssfr_prof = sfr_prof / sm_prof

    plot_profile(rplot, ssfr_prof, results_dir+'/ssfr_profile_gal_'+str(i)+'.png', 'sSFR', title=title)
    plot_profile(rplot, h1_prof, results_dir+'/h1_profile_gal_'+str(i)+'.png', 'HI mass surface density', title=title)
    plot_profile(rplot, h2_prof, results_dir+'/h2_profile_gal_'+str(i)+'.png', 'H2 mass surface density', title=title)
    plot_profile(rplot, temp_prof, results_dir+'/temp_profile_gal_'+str(i)+'.png', 'Temperature', title=title)

