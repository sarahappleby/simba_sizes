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

factor = 5.
dr = 0.2
n = int(factor / dr)
rplot = np.arange(0., dr*n, dr) + (dr*0.5)

mass_bins = [10., 10.5, 11.]
bin_labels = ['10.0 - 10.5', '10.5 - 11.0','> 11.0']
masks = [2, 3, 4]

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
if wind == 's50j7k':
    halflight_file = '/home/sapple/simba_sizes/sizes/data/halfradius_R.h5'
else:
    halflight_file = '/home/sapple/simba_sizes/sizes/data/halfradius_agn_R.h5'
snapfile = data_dir+'snap_'+model+'_'+snap+'.hdf5'

sim =  caesar.load(data_dir+'Groups/'+model+'_'+snap+'.hdf5')

h = sim.simulation.hubble_constant
redshift = sim.simulation.redshift
cosmo = FlatLambdaCDM(H0=100*h, Om0=sim.simulation.omega_matter, Ob0=sim.simulation.omega_baryon, Tcmb0=2.73)
thubble = cosmo.age(redshift).value # in Gyr

# load in the caesar galaxy data to make an initial cut of star forming galaxies
gal_cent = np.array([i.central for i in sim.galaxies])
if centrals:
    print 'Looking at centrals'
else:
    print 'Looking at satellites'
    gal_cent = np.invert(gal_cent)
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

gas_pos = readsnap(snapfile, 'pos', 'gas', suppress=1, units=1) / (h*(1.+redshift)) # in kpc
gas_mass = readsnap(snapfile, 'mass', 'gas', suppress=1, units=1) / h # in Mo
gas_temp = readsnap(snapfile, 'u', 'gas', suppress=1, units=1)

bh_pos = readsnap(snapfile, 'pos', 'bndry', suppress=1, units=1) / (h*(1.+redshift)) # in kpc

no_gals = np.zeros(len(bin_labels))

for j, m in enumerate(masks):
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
        use_gas_temp = np.zeros((len(gal_ids_use), n))

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
                For the gas particles:
                """
                pos = gas_pos[glist]
                mass = gas_mass[glist]
                temp = gas_temp[glist]

                if len(sim.galaxies[gal_ids_use[i]].bhlist) > 0.:
                        pos -= bh_pos[sim.galaxies[gal_ids_use[i]].bhlist[0]]
                else:
                        s_pos = star_pos[slist]
                        s_mass = star_mass[slist]
                        pos -= center_of_quantity(s_pos, s_mass)
                        print 'No black holes to center on, centering on stars'


                r  = np.sqrt(pos[:, 0]**2 +pos[:, 1]**2) / rhalf
                mass_profile = make_profile(n, dr, r, mass, rhalf)
                temp_profile = make_profile(n, dr, r, temp*mass, rhalf)
                use_gas_temp[i] = temp_profile / mass_profile

        with h5py.File(results_dir+'/mask_'+str(m)+'_gas_temp.h5', 'a') as f:
            f.create_dataset('temp', data=np.array(use_gas_temp))
