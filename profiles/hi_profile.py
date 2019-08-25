from scipy.ndimage.filters import gaussian_filter

import h5py
import sys
import numpy as np
import os
import gc
import matplotlib.pyplot as plt

sys.path.append('/home/sapple/tools/')
from projection import *
import plotmedian as pm

from readgadget import readsnap
import caesar
from astropy.cosmology import FlatLambdaCDM

from profile_methods import *

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=16)

dr = 0.5 # kpc
h1_limit = 1.2e6 # Msun/kpc**2
Npixels = 512
n_min = 50

model = sys.argv[1]
wind = sys.argv[2]
snap = sys.argv[3]
results_dir = sys.argv[4]

results_dir += '/'+model + '_' + snap + '/'
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

gal_h1m = np.array([i.masses['HI'].in_units('Msun') for i in sim.central_galaxies])
gal_sm = np.array([i.masses['stellar'].in_units('Msun') for i in sim.central_galaxies])
gal_sfr = np.array([i.sfr.in_units('Msun/yr') for i in sim.central_galaxies])
gal_ssfr = np.log10(gal_sfr / gal_sm)
gal_rad = np.array([i.radii['stellar_half_mass'].in_units('kpc') for i in sim.central_galaxies])
richness = np.log10(gal_h1m / gal_sm)

mass_mask = np.log10(gal_sm) > 9.
n_mask = np.array([(len(i.glist) > n_min) for i in sim.central_galaxies])
mask = mass_mask * n_mask * (gal_h1m > 0.)
mask = mass_mask * (gal_h1m > 0.)
gal_ids = np.arange(len(sim.central_galaxies))[mask]

gas_pos = readsnap(snapfile, 'pos', 'gas', suppress=1, units=1) / (h*(1.+redshift)) # in kpc
gas_mass = readsnap(snapfile, 'mass', 'gas', suppress=1, units=1) / h # in Mo
gas_h1 = readsnap(snapfile, 'NeutralHydrogenAbundance', 'gas', suppress=1, units=1)

bh_pos = readsnap(snapfile, 'pos', 'bndry', suppress=1, units=1) / (h*(1.+redshift)) # in kpc

h1_radii = np.zeros(len(gal_ids))

for i in range(len(gal_ids)):
    print('\n')
    print('Galaxy ' +str(gal_ids[i]))

    rhalf = gal_rad[gal_ids[i]]
    title = 'log M* = ' + str(round(np.log10(gal_sm[gal_ids[i]]), 2)) + '; log(sSFR) = ' + format(round(gal_ssfr[gal_ids[i]], 2)) + '; rhalf = ' + str(round(rhalf, 2)) 
       
    glist = sim.central_galaxies[gal_ids[i]].halo.glist
    pos = gas_pos[glist]
    mass = gas_mass[glist]
    h1 = gas_h1[glist]

    if len(sim.central_galaxies[gal_ids[i]].bhlist) > 0.:
        pos -= bh_pos[sim.central_galaxies[gal_ids[i]].bhlist[0]]
    else:
        pos -= center_of_quantity(pos, mass)

    r  = np.sqrt(pos[:, 0]**2 +pos[:, 1]**2) 

    h1_profile, h1_radius = hi_profile(r, dr, h1*mass, rhalf, h1_limit)
    h1_radii[i] = h1_radius[-1]

    rmask = r < rhalf*8. 
    xmin = -20.
    xmax = 20.
    im,xedges,yedges=np.histogram2d(pos[:, 0][rmask],pos[:, 1][rmask],bins=(Npixels,Npixels),weights=(h1*mass)[rmask])
    im=im/((xmax-xmin)/float(Npixels))**2
    sigma = 0.5*Npixels/(xedges.max()-xedges.min())
    im = gaussian_filter(im, sigma)
    extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
    extent = [-1.*rhalf*7., rhalf*7., -1.*rhalf*7., rhalf*7.]
    v_min = 3.
    v_max = 9.
    plt.imshow(np.log10(im.transpose()+0.0001),extent=extent, interpolation='nearest',cmap='magma',
                        vmin=v_min,vmax=v_max, origin="lower")
    plt.title(title)
    plt.colorbar(label=r'$\textrm{HI surface density}$')
    plt.savefig(results_dir+'/images/gal_'+str(gal_ids[i])+'.png')
    plt.clf()

    plt.plot(h1_radius, np.array(h1_profile)/1.e6, linestyle='--', marker='.')
    plt.axhline(1., linestyle='--', color='m')
    plt.ylabel(r'$\textrm{HI surface density}$')
    plt.xlabel(r'$R (\textrm{kpc})$')
    plt.title(title)
    plt.xlim(0, )
    plt.savefig(results_dir+'/profiles/h1_profile_gal_'+str(gal_ids[i])+'.png')
    plt.clf()


with h5py.File(results_dir+'/'+model+'_'+wind+'_'+snap+'_h1_data.h5', 'a') as f:
    f.create_dataset('gal_ids', data=np.array(gal_ids))
    f.create_dataset('h1_radii', data=np.array(h1_radii))
    f.create_dataset('h1_mass', data=np.array(gal_h1m[mask]))

bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(np.log10(gal_h1m[mask]),h1_radii)
plt.scatter(np.log10(gal_h1m[mask]), h1_radii, c=richness[mask], s=1., cmap='jet_r')
plt.colorbar(label=r'$\textrm{log} (M_{HI} / M_*)$')
plt.plot(bin_cent, ymean, '--', lw=3, color='g')
plt.xlabel(r'$\textrm{log} (M_{HI} / M_{\odot})$')
plt.ylabel(r'$R_{HI} (\textrm{kpc})$')
plt.savefig(results_dir+'/'+model+'_'+wind+'_'+snap+'_h1_size_mass.png')
plt.clf()


