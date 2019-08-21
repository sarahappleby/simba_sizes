import caesar
import sys
import h5py
import pylab as plt
import numpy as np
from readgadget import *
sys.path.append('/home/sapple/tools')
import plotmedian as pm
import pygad as pg
from astropy.cosmology import FlatLambdaCDM

plotvar = 'mstar'
plotvar = 'R'  # options are 'bolometric', 'U', 'B', 'V', 'R', and 'K'
frac = 0.5

# define input file
MODEL = 'm100n1024'
WIND = ['s50j7k', ]

snaps = ['151', '125', '105', '090', '078', '062']
snaps = ['151']
save_file = '/home/sapple/simba_sizes/sizes/data/halfradius_3d_R.h5'

def center_of_weight(pos, weight):
        if len(pos.shape) == 2:
                weight =  np.transpose(np.array([weight,]*len(pos[0])))
        return np.sum(pos*weight, axis=0) / np.sum(weight, axis=0)

def compute_3d_rfrac(mass,pos,frac=0.5):
    mtot = np.sum(mass)
    r2 = np.sqrt(np.sum(pos**2, axis=1))
    sortindex = np.argsort(r2)
    r2sort = r2[sortindex]
    msum = np.cumsum(mass[sortindex])
    for i in range(len(mass)):
        if msum[i] > frac*mtot:
            if r2sort[i] > 100*100: rh = mygals[igal].radii['stellar_half_mass']
            else: rh = np.sqrt(r2sort[i])
            break
    return rh


# load in input file
#fig,ax = plt.subplots()
for SNAP in snaps:
    for iwind in range(0,len(WIND)):
        snapfile = '/home/rad/data/'+MODEL+'/'+WIND[iwind]+'/snap_'+MODEL+'_'+SNAP+'.hdf5'
        infile = '/home/rad/data/'+MODEL+'/'+WIND[iwind]+'/Groups/'+MODEL+'_'+SNAP+'.hdf5'
        sim = caesar.load(infile,LoadHalo=False)
        redshift = np.round(sim.simulation.redshift,3)
        h = sim.simulation.hubble_constant

        mygals = sim.galaxies
        print 'Doing z=',redshift,'snapshot',snapfile,'with',len(mygals),'galaxies.'

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
        bh_pos = readsnap(snapfile, 'pos', 'bndry', suppress=1, units=1) / (h*(1.+redshift)) # in kpc

        rhalf = np.zeros(len(mygals))
        rhalfmass = np.zeros(len(mygals))

        
        pgsnap = pg.Snap(snapfile,physical=True)
        print 'Computing',plotvar,'magnitudes for',len(pgsnap.stars),'stars.'
        allL = pg.snapshot.get_luminosities(pgsnap.stars,band=plotvar)
        
        for igal in range(len(mygals)):
                mass = np.array([sm[k] for k in mygals[igal].slist])
                lum = np.array([allL[k] for k in mygals[igal].slist])
                pos = np.array([sp[k] for k in mygals[igal].slist])
                if len(mygals[igal].bhlist) > 0.:
                    pos -= bh_pos[mygals[igal].bhlist[0]]
                else:
                    pos -= center_of_weight(pos, mass)
                rhalfmass[igal] = compute_3d_rfrac(mass,pos, frac)
                rhalf[igal] = compute_3d_rfrac(lum,pos, frac)

                if igal%(len(mygals)/20)==0: print('%d logM*= %.3f rh2d_new= %.2f rh2dMs_new= %.2f rh3d= %.3f'%(igal,np.log10(ms[igal]),rhalf[igal],rhalfmass[igal],mygals[igal].radii['stellar_half_mass']))
        
        with h5py.File(save_file, 'a') as hf:
            hf.create_dataset(MODEL+'_'+WIND[iwind]+'_'+str(SNAP)+'_halflight', data=np.array(rhalf))
            hf.create_dataset(MODEL+'_'+WIND[iwind]+'_'+str(SNAP)+'_halfmass', data=np.array(rhalfmass))


