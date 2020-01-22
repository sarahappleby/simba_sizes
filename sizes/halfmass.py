
import caesar
import sys
import h5py
import numpy as np
from readgadget import *
import pygad as pg
from astropy.cosmology import FlatLambdaCDM

from size_finding_methods import compute_rfrac

plotvar = 'mstar'

# define input file
MODEL = 'm25n512'
WIND = ['s50j7k', 's50noagn', 's50nojet', 's50nox']
WIND = ['s50']
mlim = 8.7
mmax = 12.7
frac = 0.5

snaps = ['151', '125', '105', '090', '078', '062']
save_file = '/home/sapple/simba_sizes/sizes/data/halfmass.h5'


for SNAP in snaps:
    for iwind in range(0,len(WIND)):
        
        snapfile = '/home/rad/data/'+MODEL+'/'+WIND[iwind]+'/snap_'+MODEL+'_'+SNAP+'.hdf5'
        infile = '/home/rad/data/'+MODEL+'/'+WIND[iwind]+'/Groups/'+MODEL+'_'+SNAP+'.hdf5'
        sim = caesar.load(infile,LoadHalo=False)
        redshift = np.round(sim.simulation.redshift,3)
        h = sim.simulation.hubble_constant
        
        mygals = sim.galaxies
        print('Doing z=',redshift,'snapshot',snapfile,'with',len(mygals),'galaxies.')
        
        # read in galaxy info
        ids = np.asarray([i.GroupID for i in mygals])
        central = np.asarray([i.central for i in mygals])
        ms = np.asarray([i.masses['stellar'] for i in mygals])
        mbh = np.asarray([i.masses['bh'] for i in mygals])
        sfr = np.asarray([i.sfr for i in mygals])

	# read in particle info
        sm = readsnap(snapfile,'mass','star',units=1,suppress=1)/h  # Mo
        sp = readsnap(snapfile,'pos','star',units=1,suppress=1)/(1+redshift)/h # pkpc
        gm = readsnap(snapfile,'mass','gas',units=1,suppress=1)/h # Mo
        gp = readsnap(snapfile,'pos','gas',units=1,suppress=1)/(1+redshift)/h # pkpc
        gsfr = readsnap(snapfile,'sfr','gas',units=1,suppress=1) # Mo/yr
        
        # compute rhalf
        rhalf = np.zeros([3,len(mygals)])
        
        for igal in range(len(mygals)):
            mass = np.array([sm[k] for k in mygals[igal].slist])
            pos = np.array([sp[k] for k in mygals[igal].slist])
            for idir0 in range(3): 
                rhalf[idir0][igal],cent = compute_rfrac(idir0,mass,pos, frac)
            if igal%20==0: print('%d logM*= %.3f %.3f c= [%.3f,%.3f] rh2d= %.3f %.3f %.3f rh3d= %.3f'%(igal,np.log10(ms[igal]),np.log10(sum(mass)),cent[0],cent[1],rhalf[0][igal],rhalf[1][igal],rhalf[2][igal],mygals[igal].radii['stellar_half_mass']))        
        with h5py.File(save_file, 'a') as f:
            f.create_dataset(MODEL+'_'+WIND[iwind]+'_'+str(SNAP), data=np.array(rhalf))

