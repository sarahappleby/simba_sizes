
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

rotate_galaxies = True
vec = np.array([0, 0, 1]) # face-on projection to collapse the z direction

factor = 5.
dr = 0.2
n = int(factor / dr)
rplot = np.arange(0., dr*n, dr) + (dr*0.5)

model = sys.argv[1]
wind = sys.argv[2]
snap = sys.argv[3]
results_dir = sys.argv[4]

results_dir += '/'+wind+'/'+model+'_'+snap+'/'
if not os.path.exists(results_dir):
        os.makedirs(results_dir)

if rotate_galaxies:
        results_file = results_dir+'all_profiles_rotated_faceon.h5'
else:
        results_file = results_dir+'all_profiles_random_orientation.h5'

data_dir = '/home/rad/data/'+model+'/'+wind+'/'
if snap == '151':
    halflight_file = '/home/sapple/simba_sizes/sizes/data/pyloser_sdss_r.h5'
else:
    halflight_file = '/home/sapple/simba_sizes/sizes/data/pyloser_v.h5'
snapfile = data_dir+'snap_'+model+'_'+snap+'.hdf5'

sim =  caesar.load(data_dir+'Groups/'+model+'_'+snap+'.hdf5')

h = sim.simulation.hubble_constant
redshift = sim.simulation.redshift
cosmo = FlatLambdaCDM(H0=100*h, Om0=sim.simulation.omega_matter, Ob0=sim.simulation.omega_baryon, Tcmb0=2.73)
thubble = cosmo.age(redshift).value # in Gyr

# load in the caesar galaxy data to make an initial cut of star forming galaxies
gal_cent = np.array([i.central for i in sim.galaxies])
gal_sm = np.array([i.masses['stellar'].in_units('Msun') for i in sim.galaxies])
gal_sfr = np.array([i.sfr.in_units('Msun/yr') for i in sim.galaxies])
gal_ssfr = gal_sfr / gal_sm
gal_sm = np.log10(gal_sm)
gal_ssfr = np.log10(gal_ssfr)

with h5py.File(halflight_file, 'r') as f:
        gal_rad = f['abs_'+model+'_'+wind+'_'+snap][:] # these are in pkpc
gal_rad = np.sum(gal_rad, axis=0) / 3.

gal_ids = np.array([True for i in range(len(sim.galaxies))])
gal_ids = np.arange(len(gal_ids))[gal_ids]

# load in star particle data with readsnap
star_pos = readsnap(snapfile, 'pos', 'star', suppress=1, units=1) / (h*(1.+redshift)) # in kpc
star_mass = readsnap(snapfile, 'mass', 'star', suppress=1, units=1) / h # in Mo
star_vels = readsnap(snapfile, 'vel', 'star', suppress=1, units=0) # in km/s
star_tform = readsnap(snapfile, 'age', 'star', suppress=1, units=1) # expansion times at time of formation

gas_pos = readsnap(snapfile, 'pos', 'gas', suppress=1, units=1) / (h*(1.+redshift)) # in kpc
gas_mass = readsnap(snapfile, 'mass', 'gas', suppress=1, units=1) / h # in Mo
gas_sfr = readsnap(snapfile, 'sfr', 'gas', suppress=1, units=1) # in Mo/yr
gas_h1 = readsnap(snapfile, 'NeutralHydrogenAbundance', 'gas', suppress=1, units=1)
gas_h2 = readsnap(snapfile, 'fh2', 'gas', suppress=1, units=1)
gas_temp = readsnap(snapfile, 'u', 'gas', suppress=1, units=1)

bh_pos = readsnap(snapfile, 'pos', 'bndry', suppress=1, units=1) / (h*(1.+redshift)) # in kpc

all_star_m = np.zeros((len(sim.galaxies), n))
all_gas_h1 = np.zeros((len(sim.galaxies), n))
all_gas_h2 = np.zeros((len(sim.galaxies), n))
all_gas_h1_m = np.zeros((len(sim.galaxies), n))
all_gas_h2_m = np.zeros((len(sim.galaxies), n))
all_gas_m = np.zeros((len(sim.galaxies), n))
all_gas_sfr = np.zeros((len(sim.galaxies), n))
all_gas_temp = np.zeros((len(sim.galaxies), n))
all_gas_npart = np.zeros((len(sim.galaxies), n))
all_star_npart = np.zeros((len(sim.galaxies), n))

for i in gal_ids:
    
    print('\n')
    print('Galaxy ' +str(i))
    
    if (sim.galaxies[i].halo is None):
        print('Galaxy has no halo, skipping')
        continue

    slist = sim.galaxies[i].halo.slist
    glist = sim.galaxies[i].halo.glist
    rhalf = gal_rad[i]
    title = 'log M* = ' + str(round(gal_sm[i], 2)) + '; log(sSFR) = ' + format(round(gal_ssfr[i], 2))

    print(str(len(glist)) + ' gas particles')
    print('log sSFR: ' + format(round(gal_ssfr[i], 2)))

    """
    Get star particle data and correct for bh center or star com
    """
    pos = star_pos[slist]
    vel = star_vels[slist]
    mass = star_mass[slist]
    
    if len(sim.galaxies[i].bhlist) > 0.:
        pos -= bh_pos[sim.galaxies[i].bhlist[0]]
    else:
        pos -= center_of_quantity(pos, mass)
        print('No black holes to center on, centering on stars')
    vel -= center_of_quantity(vel, mass)
    
    if rotate_galaxies:
        axis, angle = compute_rotation_to_vec(pos[:, 0], pos[:, 1], pos[:, 2], vel[:, 0], vel[:, 1], vel[:, 2], mass, vec)
        pos[:, 0], pos[:, 1], pos[:, 2] = rotate(pos[:, 0], pos[:, 1], pos[:, 2], axis, angle)
        
    r = np.sqrt(pos[:, 0]**2 +pos[:, 1]**2) / rhalf

    all_star_m[i] = make_profile(n, dr, r, mass, rhalf)
    all_star_npart[i] = npart_profile(n, dr, r)
    
    del pos, vel, mass, r

    if len(glist) > 0.:
        """
        For the gas particles:
        """
        pos = gas_pos[glist]
        mass = gas_mass[glist]
        sfr = gas_sfr[glist]
        h1 = gas_h1[glist]
        h2 = gas_h2[glist]
        temp = gas_temp[glist]
        
        if len(sim.galaxies[i].bhlist) > 0.:
            pos -= bh_pos[sim.galaxies[i].bhlist[0]]
        else:
            s_pos = star_pos[slist]
            s_mass = star_mass[slist]
            pos -= center_of_quantity(s_pos, s_mass)
            print('No black holes to center on, centering on stars')
            
        if rotate_galaxies:
            pos[:, 0], pos[:, 1], pos[:, 2] = rotate(pos[:, 0], pos[:, 1], pos[:, 2], axis, angle)
        
        r  = np.sqrt(pos[:, 0]**2 +pos[:, 1]**2) / rhalf
        
        all_gas_m[i] = make_profile(n, dr, r, mass, rhalf)
        all_gas_sfr[i] = make_profile(n, dr, r, sfr, rhalf)
        all_gas_h1[i] = make_profile(n, dr, r, h1, rhalf)
        all_gas_h2[i] = make_profile(n, dr, r, h2, rhalf)
        all_gas_temp[i] = make_profile(n, dr, r, temp*mass, rhalf)
        all_gas_temp[i] /= all_gas_m[i]
        all_gas_h1_m[i] = make_profile(n, dr, r, h1*mass, rhalf)
        all_gas_h2_m[i] = make_profile(n, dr, r, h2*mass, rhalf)
        all_gas_npart[i] = npart_profile(n, dr, r) 

        del pos, mass, sfr, h1, h2, temp, r

    gc.collect()


with h5py.File(results_file, 'a') as f:
    f.create_dataset('gas_sfr', data=np.array(all_gas_sfr))
    f.create_dataset('h1', data=np.array(all_gas_h1))
    f.create_dataset('h2', data=np.array(all_gas_h2))
    f.create_dataset('gm', data=np.array(all_gas_m))
    f.create_dataset('sm', data=np.array(all_star_m))
    f.create_dataset('temp', data=np.array(all_gas_temp))
    f.create_dataset('h1_m', data=np.array(all_gas_h1_m))
    f.create_dataset('h2_m', data=np.array(all_gas_h2_m))
    f.create_dataset('gas_npart', data=np.array(all_gas_npart))
    f.create_dataset('star_npart', data=np.array(all_star_npart))
    f.create_dataset('gal_ids', data=np.array(gal_ids))
