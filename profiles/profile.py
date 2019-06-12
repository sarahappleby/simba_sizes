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

from profile_methods import *

age_min = 0.
age_max = 150.
mass_loss = 1.18
time = (age_max - age_min) * 1.e6
rotate_galaxies = False
vec = np.array([0, 0, 1]) # face-on projection to collapse the z direction

factor = 5.
dr = 0.2
n = int(factor / dr)
rplot = np.arange(0., dr*n, dr) + (dr*0.5)

# for making profiles
mass_bins = [9., 9.5, 10., 10.5, 11.]
bin_labels = ['9.0 - 9.5', '9.5 - 10.0', '10.0 - 10.5', '10.5 - 11.0', '> 11.0']
masks = range(len(bin_labels))

model = sys.argv[1]
wind = sys.argv[2]
snap = sys.argv[3]
sample_file = sys.argv[4]
results_dir = sys.argv[5]

if 'gv' in sample_file.split('/', -1)[-1]:
    selection = 'green_valley'
elif 'sf' in sample_file.split('/', -1)[-1]:
    selection = 'star_forming'

if 'satellites' in results_dir:
    centrals = False
    print 'Looking at satellite galaxies'
else:
    centrals = True
    print 'Looking at central galaxies'

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

# load in the caesar galaxy data to make an initial cut of star forming galaxies
gal_cent = np.array([i.central for i in sim.galaxies])
if centrals: gal_cent = np.invert(gal_cent)
gal_sm = np.array([i.masses['stellar'].in_units('Msun') for i in sim.galaxies])
gal_sfr = np.array([i.sfr.in_units('Msun/yr') for i in sim.galaxies])
gal_ssfr = gal_sfr / gal_sm

with h5py.File(halflight_file, 'r') as f:
    gal_rad = f[model+'_'+wind+'_'+snap+'_halflight'][:] # these are in pkpc
gal_rad = np.sum(gal_rad, axis=0) / 3.
#gal_rad = np.array([i.radii['stellar_half_mass'].in_units('kpc') for i in sim.galaxies])

gal_sm = np.log10(gal_sm[gal_ids*gal_cent])
gal_ssfr = np.log10(gal_ssfr[gal_ids*gal_cent])
gal_rad = gal_rad[gal_ids*gal_cent]
gal_ids = np.arange(len(gal_ids))[gal_ids*gal_cent]

# load in star particle data with readsnap
star_pos = readsnap(snapfile, 'pos', 'star', suppress=1, units=1) / (h*(1.+redshift)) # in kpc
star_mass = readsnap(snapfile, 'mass', 'star', suppress=1, units=1) / h # in Mo
star_vels = readsnap(snapfile, 'vel', 'star', suppress=1, units=0) # in km/s

gas_pos = readsnap(snapfile, 'pos', 'gas', suppress=1, units=1) / (h*(1.+redshift)) # in kpc
gas_mass = readsnap(snapfile, 'mass', 'gas', suppress=1, units=1) / h # in Mo
gas_vels = readsnap(snapfile, 'vel', 'gas', suppress=1, units=0) # in km/s
gas_sfr = readsnap(snapfile, 'sfr', 'gas', suppress=1, units=1) # in Mo/yr
gas_h1 = readsnap(snapfile, 'NeutralHydrogenAbundance', 'gas', suppress=1, units=1)
gas_h2 = readsnap(snapfile, 'fh2', 'gas', suppress=1, units=1)

bh_pos = readsnap(snapfile, 'pos', 'bndry', suppress=1, units=1) / (h*(1.+redshift)) # in kpc

no_gals = np.zeros(len(bin_labels))

for m in masks:
        print '\n'
        print 'Looking at mass bin ' + bin_labels[m]
        if m != len(mass_bins) - 1:
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

        use_star_m = np.zeros((len(gal_ids_use), n))
        use_gas_m = np.zeros((len(gal_ids_use), n))
        use_gas_sfr = np.zeros((len(gal_ids_use), n))
        use_gas_h1 = np.zeros((len(gal_ids_use), n))
        use_gas_h2 = np.zeros((len(gal_ids_use), n))


        for i in range(len(gal_ids_use)):

                print '\n'
                print 'Galaxy ' +str(gal_ids_use[i])
                slist = sim.galaxies[gal_ids_use[i]].halo.slist
                glist = sim.galaxies[gal_ids_use[i]].halo.glist
                rhalf = gal_rad_use[i]
                title = 'log M* = ' + str(round(gal_sm_use[i], 2)) + '; log(sSFR) = ' + format(round(gal_ssfr_use[i], 2))

                print str(len(glist)) + ' gas particles'
                print 'log sSFR: ' + format(round(gal_ssfr_use[i], 2))

                """
                Get star particle data and correct for bh center or star com
                """
                pos = star_pos[slist]
                vel = star_vels[slist]
                mass = star_mass[slist]
                if len(sim.galaxies[gal_ids_use[i]].bhlist) > 0.:
                        pos -= bh_pos[sim.galaxies[gal_ids_use[i]].bhlist[0]]
                else:
                        pos -= center_of_quantity(pos, mass)
                        print 'No black holes to center on, centering on stars'
                vel -= center_of_quantity(vel, mass)

                if rotate_galaxies:
                        axis, angle = compute_rotation_to_vec(pos[:, 0], pos[:, 1], pos[:, 2], vel[:, 0], vel[:, 1], vel[:, 2], mass, vec)
                        pos[:, 0], pos[:, 1], pos[:, 2] = rotate(pos[:, 0], pos[:, 1], pos[:, 2], axis, angle)

                r = np.sqrt(pos[:, 0]**2 +pos[:, 1]**2) / rhalf

                plot_name = results_dir + '/images/gal_'+str(gal_ids_use[i]) + '_stars.png'
                make_image(pos[:, 0], pos[:, 1], np.log10(mass), plot_name, rhalf)

                use_star_m[i] = make_profile(n, dr, r, mass, rhalf)

                """
                For the gas particles:
                """
                pos = gas_pos[glist]
                vel = gas_vels[glist]
                mass = gas_mass[glist]
                sfr = gas_sfr[glist]
                h1 = gas_h1[glist]
                h2 = gas_h2[glist]

                if len(sim.galaxies[gal_ids_use[i]].bhlist) > 0.:
                        pos -= bh_pos[sim.galaxies[gal_ids_use[i]].bhlist[0]]
                else:
                        pos -= center_of_quantity(pos, mass)
                        print 'No black holes to center on, centering on stars'
                vel -= center_of_quantity(vel, mass)

                if rotate_galaxies:
                        pos[:, 0], pos[:, 1], pos[:, 2] = rotate(pos[:, 0], pos[:, 1], pos[:, 2], axis, angle)

                r  = np.sqrt(pos[:, 0]**2 +pos[:, 1]**2) / rhalf

                plot_name = results_dir + '/images/gal_'+str(gal_ids_use[i]) + '_gas.png'
                make_image(pos[:, 0], pos[:, 1], np.log10(mass), plot_name, rhalf)

                use_gas_m[i] = make_profile(n, dr, r, mass, rhalf)
                use_gas_sfr[i] = make_profile(n, dr, r, sfr, rhalf)
                use_gas_h1[i] = make_profile(n, dr, r, h1, rhalf)
                use_gas_h2[i] = make_profile(n, dr, r, h2, rhalf)
                
                plot_profile(rplot, use_gas_sfr[i], results_dir+'/profiles/gas_sfr_profile_gal_'+str(gal_ids_use[i])+'.png', 'SFR surface density', title=title)
                plot_profile(rplot, use_gas_h1[i], results_dir+'/profiles/gas_h1_profile_gal_'+str(gal_ids_use[i])+'.png', 'HI fraction surface density', title=title)
                plot_profile(rplot, use_gas_h2[i], results_dir+'/profiles/gas_h2_profile_gal_'+str(gal_ids_use[i])+'.png', 'HII fraction surface density', title=title)

        with h5py.File(results_dir+'/mask_'+str(m)+'_all_profiles.h5', 'a') as f:
                f.create_dataset('gas_sfr', data=np.array(use_gas_sfr))
                f.create_dataset('h1', data=np.array(use_gas_h1))
                f.create_dataset('h2', data=np.array(use_gas_h2))
                f.create_dataset('gm', data=np.array(use_gas_m))
                f.create_dataset('sm', data=np.array(use_star_m))
                f.create_dataset('gal_ids', data=np.array(gal_ids_use))

        del use_star_m, use_gas_sfr, use_gas_h1, use_gas_h2
        gc.collect()
        print '\n'

