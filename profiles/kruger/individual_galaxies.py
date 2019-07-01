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

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 14})

dr = 1. # kpc
factor = 5.
ssfr_limit = -1.8
xlabel = r'$R (kpc)$'
clabel = r'$ \textrm{log} (\Sigma_{*} / M_{\odot}\textrm{kpc}^{-2})$'
mass_bins = [10., 10.5, 11.]
bin_labels = ['10.0 - 10.5', '10.5 - 11.0']

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
cosmo = FlatLambdaCDM(H0=100*h, Om0=sim.simulation.omega_matter, Ob0=sim.simulation.omega_baryon, Tcmb0=2.73)
thubble = cosmo.age(redshift).value # in Gyr

gal_n = np.array([len(i.glist) for i in sim.central_galaxies])
gal_ids = np.arange(len(sim.central_galaxies))
gal_rad = np.array([i.radii['stellar_half_mass'].in_units('kpc') for i in sim.central_galaxies])
gal_sm = np.array([i.masses['stellar'].in_units('Msun') for i in sim.central_galaxies])
gal_sfr = np.array([i.sfr.in_units('Msun/Gyr') for i in sim.central_galaxies])
gal_ssfr = np.log10(gal_sfr / gal_sm)
gal_sm = np.log10(gal_sm)

sf_mask = gal_ssfr > ssfr_limit
n_mask = gal_n > 1000

star_pos = readsnap(snapfile, 'pos', 'star', suppress=1, units=1) / (h*(1.+redshift)) # in kpc
star_mass = readsnap(snapfile, 'mass', 'star', suppress=1, units=1) / h # in Mo
star_tform = readsnap(snapfile, 'age', 'star', suppress=1, units=1) # expansion times at time of formation

gas_pos = readsnap(snapfile, 'pos', 'gas', suppress=1, units=1) / (h*(1.+redshift)) # in kpc
gas_mass = readsnap(snapfile, 'mass', 'gas', suppress=1, units=1) / h # in Mo
gas_sfr = readsnap(snapfile, 'sfr', 'gas', suppress=1, units=1) # in Mo/yr
gas_h1 = readsnap(snapfile, 'NeutralHydrogenAbundance', 'gas', suppress=1, units=1)
gas_h2 = readsnap(snapfile, 'fh2', 'gas', suppress=1, units=1)
gas_temp = readsnap(snapfile, 'u', 'gas', suppress=1, units=1)

bh_pos = readsnap(snapfile, 'pos', 'bndry', suppress=1, units=1) / (h*(1.+redshift)) # in kpc

gals = [45, 94, 195, 201, 259, 364]
"""
gals = []

for j, b in enumerate(bin_labels):
    print '\n'
    print 'Looking at mass bin ' + b

    if j != len(mass_bins) - 1:
        sm_mask = (gal_sm > mass_bins[j]) & (gal_sm < mass_bins[j+1])
    else:
        sm_mask = gal_sm > mass_bins[j]

    # star forming:

    gal_ids_use = gal_ids[sm_mask*sf_mask*n_mask]
    if len(gal_ids_use) > 0.:
        choose = np.random.randint(len(gal_ids_use))
        gals.append(gal_ids_use[choose])
        print 'Found star forming galaxy'

    # quenched galaxies:

    gal_ids_use = gal_ids[sm_mask*np.invert(sf_mask)*n_mask]
    if len(gal_ids_use) > 0.:
        choose = np.random.randint(len(gal_ids_use))
        gals.append(gal_ids_use[choose])
        print 'Found quenched galaxy'
"""
for i in gals:
    print 'Galaxy ' +str(i)
    slist = sim.central_galaxies[i].halo.slist
    glist = sim.central_galaxies[i].halo.glist
    rhalf = gal_rad[i]
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
    age_prof = real_profile(n, dr, r, ages*mass)
    age_prof /= sm_prof
    plot_profile(rplot, np.log10(sm_prof), results_dir+'/sm_profile_gal_'+str(i)+'.png', r'$ \textrm{log} (\Sigma_{*} / M_{\odot}\textrm{kpc}^{-2})$', xlabel=xlabel, title=title)
    plot_profile(rplot, np.log10(age_prof), results_dir+'/age_profile_gal_'+str(i)+'.png', r'$ \textrm{log} (\tau_{*} / \textrm{Myr})$', xlabel=xlabel, title=title)

    plot_name = results_dir + '/gal_'+str(i) + '_stars.png'
    make_image(pos[r < rhalf*factor][:, 0], pos[r < rhalf*factor][:, 1], np.log10(mass)[r < rhalf*factor], plot_name, rhalf, clabel=clabel)

    pos = gas_pos[glist]
    mass = gas_mass[glist]
    sfr = gas_sfr[glist]
    h1 = gas_h1[glist]
    h2 = gas_h2[glist]
    temp = gas_temp[glist]

    if len(sim.central_galaxies[i].bhlist) > 0.:
        pos -= bh_pos[sim.central_galaxies[i].bhlist[0]]
    else:
        pos -= center_of_quantity(pos, mass)
        print 'No black holes to center on, centering on stars'

    r = np.sqrt(pos[:, 0]**2 +pos[:, 1]**2)

    gas_m_prof = real_profile(n, dr, r, mass)
    sfr_prof = real_profile(n, dr, r, sfr)
    h1_prof = real_profile(n, dr, r, h1*mass)
    h2_prof = real_profile(n, dr, r, h2*mass)
    temp_prof = real_profile(n, dr, r, temp*mass)
    temp_prof /= gas_m_prof
    ssfr_prof = sfr_prof / sm_prof

    plot_profile(rplot, np.log10(ssfr_prof), results_dir+'/ssfr_profile_gal_'+str(i)+'.png', r'$\textrm{log} (\textrm{sSFR} / \textrm{Gyr}^{-1})$', xlabel=xlabel, title=title)
    plot_h1_profile(rplot, h1_prof / 1.e6, results_dir+'/h1_profile_gal_'+str(i)+'.png', r'$ \Sigma_{HI} (M_{\odot}\textrm{pc}^{-2})$', xlabel=xlabel, title=title)
    plot_profile(rplot, h2_prof / 1.e6, results_dir+'/h2_profile_gal_'+str(i)+'.png', r'$\Sigma_{H_2} (M_{\odot}\textrm{pc}^{-2})$', xlabel=xlabel, title=title)
    plot_profile(rplot, np.log10(temp_prof), results_dir+'/temp_profile_gal_'+str(i)+'.png', r'$\textrm{log} (T /K)$', xlabel=xlabel, title=title)

    plot_name = results_dir + '/gal_'+str(i) + '_h1.png'
    make_image(pos[r < rhalf*factor][:, 0], pos[r < rhalf*factor][:, 1], np.log10(mass*h1)[r < rhalf*factor], plot_name, rhalf, clabel=clabel)
    #plot_name = results_dir + '/gal_'+str(i) + '_h2.png'
    #make_image(pos[:, 0][r < rhalf*factor], pos[:, 1][r < rhalf*factor], np.log10(mass*h2)[r < rhalf*factor], plot_name, rhalf, clabel=clabel)


