
import caesar
import sys
import h5py
import pylab as plt
import numpy as np
import pygad as pg
from astropy.cosmology import FlatLambdaCDM
from readgadget import readsnap
sys.path.append('/home/sapple/tools')
import plotmedian as pm

from size_finding_methods import *

if __name__ == '__main__':
    
    mtype = 'abs'
    band = 'v'
    band = 'sdss_r'
    save_file = '/home/sapple/simba_sizes/sizes/data/pyloser_'+band+'.h5'

    # define input file
    MODEL = 'm50n512'
    WIND = ['s50j7k', 's50noagn', 's50nojet', 's50nox']
    WIND = ['s50noagn']
    mlim = 8.7
    mmax = 12.7
    snaps = ['151',]

    for SNAP in snaps:
        
        for iwind in range(0,len(WIND)):
            
            loserfile = '/disk01/rad/sim/'+MODEL+'/'+WIND[iwind]+'/Groups/pyloser_stars_'+MODEL+'_'+SNAP+'.hdf5'

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
            if (MODEL == 'm100n1024') & (SNAP in ['105', '151']):
                params, mags = get_mags_band(loserfile, mtype, band, only_mtype=True, verbose=True)
                where = int(np.where(params['bands'] == band)[0])
                mags = mags.transpose()[where]
            else:
                params, mags = get_mags_band(loserfile, mtype, band, verbose=True)
            allL = absmag_to_lum(mags)
            
            slist_start = 0 
            for igal in range(len(mygals)):
                slist = mygals[igal].slist
                slist_use = list(range(slist_start, slist_start + len(slist)))
                lum = np.array(allL[slist_use])
                pos = np.array(sp[slist])
                slist_start += len(slist)
                for idir0 in range(3): 
                    rhalf[idir0][igal],cent = compute_rfrac(idir0,lum,pos)
                if igal%20==0: 
                    print('%d logM*= %.3f c= [%.3f,%.3f] rh2d= %.2f %.2f %.2f rh3d= %.3f'
                            %(igal,np.log10(ms[igal]),cent[0],cent[1],rhalf[0][igal],rhalf[1][igal],rhalf[2][igal],mygals[igal].radii['stellar_half_mass']))
                    
            with h5py.File(save_file, 'a') as hf:
                hf.create_dataset(mtype+'_'+MODEL+'_'+WIND[iwind]+'_'+str(SNAP), data=np.array(rhalf))
