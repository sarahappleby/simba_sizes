import h5py
import sys
import numpy as np
import os
import gc
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM

sys.path.append('/home/sapple/tools/')
from projection import *

from readgadget import readsnap
import caesar

from profile_methods import *

# some constants and running options:
rotate_galaxies = False
factor = 3.
dr = 0.2
n = int(factor / dr)
rplot = np.arange(0., dr*n, dr) + (dr*0.5)
mass_bins = [10., 10.5, 11.]
bin_labels = ['10.0 - 10.5', '10.5 - 11.0', '> 11.0']
masks = [2, 3, 4]

# input and housekeeping:
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
else:
    centrals = True

results_dir += '/'+model + '_' + snap + '/' + wind + '/' + selection
if rotate_galaxies:
        results_dir += '/rotated_faceon'
else:
        results_dir += '/random_orientation'
if not os.path.exists(results_dir):
        os.makedirs(results_dir)
        os.makedirs(results_dir+'/images')
        os.makedirs(results_dir+'/profiles')

# get ids of galaxies to use:
with h5py.File(sample_file, 'r') as f:
        try:
                gal_ids = f[model+'_'+snap].value
        except KeyError:
                print 'Need to identify galaxies; run gv_sample.py first'

# load in caesar file:
data_dir = '/home/rad/data/'+model+'/'+wind+'/'
halflight_file = '/home/sapple/simba_sizes/sizes/data/halfradius_R.h5'
snapfile = data_dir+'snap_'+model+'_'+snap+'.hdf5'
sim =  caesar.load(data_dir+'Groups/'+model+'_'+snap+'.hdf5', LoadHalo=False)

# get cosmology info:
h = sim.simulation.hubble_constant
redshift = sim.simulation.redshift
cosmo = FlatLambdaCDM(H0=100*h, Om0=sim.simulation.omega_matter, Ob0=sim.simulation.omega_baryon, Tcmb0=2.73)
thubble = cosmo.age(redshift).value # in Gyr

# load in the caesar galaxy data:
gal_cent = np.array([i.central for i in sim.galaxies])
if centrals:
    print 'Looking at centrals'
else:
    print 'Looking at satellites'
    gal_cent = np.invert(gal_cent)
gal_sm = np.array([i.masses['stellar'].in_units('Msun') for i in sim.galaxies])
gal_sfr = np.array([i.sfr.in_units('Msun/yr') for i in sim.galaxies])
gal_ssfr = gal_sfr / gal_sm

# load in half light radii:
with h5py.File(halflight_file, 'r') as f:
    gal_rad = f[model+'_'+wind+'_'+snap+'_halflight'][:] # these are in pkpc
gal_rad = np.sum(gal_rad, axis=0) / 3.

# select central galaxies:
gal_sm = np.log10(gal_sm[gal_ids*gal_cent])
gal_ssfr = np.log10(gal_ssfr[gal_ids*gal_cent])
gal_rad = gal_rad[gal_ids*gal_cent]
gal_ids = np.arange(len(gal_ids))[gal_ids*gal_cent]

# load in star particle data with readsnap
star_pos = readsnap(snapfile, 'pos', 'star', suppress=1, units=1) / (h*(1.+redshift)) # in kpc
star_tform = readsnap(snapfile, 'age', 'star', suppress=1, units=1) # expansion times at time of formation
star_mass = readsnap(snapfile, 'mass', 'star', suppress=1, units=1) / h # in Mo
bh_pos = readsnap(snapfile, 'pos', 'bndry', suppress=1, units=1) / (h*(1.+redshift)) # in kpc

no_gals = np.zeros(len(bin_labels))

for j, m in enumerate(masks):
        
        # for each mass bin make a mass mask and select those galaxies:

        print '\n'
        print 'Looking at mass bin ' + bin_labels[j]
        if j != len(mass_bins) - 1:
                sm_mask = (gal_sm > mass_bins[j]) & (gal_sm < mass_bins[j+1])
        else:
                sm_mask = gal_sm > mass_bins[j]

        gal_sm_use = gal_sm[sm_mask]
        gal_rad_use = gal_rad[sm_mask]
        gal_ssfr_use = gal_ssfr[sm_mask]
        gal_ids_use = gal_ids[sm_mask]

        no_gals[j] = len(gal_ids_use)
        print str(no_gals[j]) + ' galaxies in bin'
        print '\n'

        use_star_ages = np.zeros((len(gal_ids_use), n))

        for i in range(len(gal_ids_use)):

                print '\n'
                print 'Galaxy ' +str(gal_ids_use[i])
        
                # get particles:
                slist = sim.galaxies[gal_ids_use[i]].slist
                rhalf = gal_rad_use[i]
               
                """
                Get star particle data and correct for bh center or star com
                """
                pos = star_pos[slist]
                mass = star_mass[slist]
                
                # compute ages:
                tform = star_tform[slist]
                print 'Getting star ages'
                ages = tage(cosmo, thubble, tform) * 1000. # get star ages in Myr

                # correct position for bh center or star com:
                if len(sim.galaxies[gal_ids_use[i]].bhlist) > 0.:
                        pos -= bh_pos[sim.galaxies[gal_ids_use[i]].bhlist[0]]
                else:
                        pos -= center_of_quantity(pos, mass)
                        print 'No black holes to center on, centering on stars'
                r = np.sqrt(pos[:, 0]**2 +pos[:, 1]**2) / rhalf
                
                # compute stellar mass weighted age profile:
                mass_profile = make_profile(n, dr, r, mass, rhalf)
                age_profile = make_profile(n, dr, r, ages*mass, rhalf)
                use_star_ages[i] = age_profile/mass_profile

        with h5py.File(results_dir+'/mask_'+str(m)+'_star_ages.h5', 'a') as f:
            f.create_dataset('ages', data=np.array(use_star_ages))

