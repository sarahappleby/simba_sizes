
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

"""
if len(sys.argv) < 3:
    print 'usage: MODEL WIND1 WIND2 ...'
    exit()
"""
plotvar = 'mstar'
plotvar = 'R'  # options are 'bolometric', 'U', 'B', 'V', 'R', and 'K'
if len(sys.argv)==5: plotvar=sys.argv[4]

# define input file
MODEL = 'm50n512'
#SNAP = int(sys.argv[2])
WIND = ['s50j7k', 's50noagn', 's50nojet', 's50nox']
WIND = ['s50nojet', 's50nox']
#plt.rc('text', usetex=True)
mlim = 8.7
mmax = 12.7
conc = False

snaps = ['151', '125', '105', '090', '078', '062']
snaps = ['151']
save_file = '/home/sapple/simba_sizes/sizes/data/halfradius_agn_R.h5'
#plot_dir = '/home/sapple/simba_sizes/sizes/plots/'

'''
if plotvar not in ['mstar','sfr']:
    try: sdss_bands = fsps.find_filter(plotvar)
    except: sys.exit('Filter %s not found'%plotvar)
    print('Doing rhalf in band %s, generating FSPS stellar pop...'%plotvar)
    spop = fsps.StellarPopulation(zcontinuous=1, sfh=0, logzsol=0.0, dust_type=2, dust2=0.2)
    #print spop.ssp_ages
'''

def compute_rfrac(idir0,mass,pos, frac=0.5):
    idir1 = (idir0+1)%3
    idir2 = (idir0+2)%3
    mtot = sum(mass)
    #cent = [sum([(mass[i]*pos[i][idir1]) for i in range(len(mass))])/mtot,sum([(mass[i]*pos[i][idir2]) for i in range(len(mass))])/mtot]  # centroid 
    cent = [sum([(pos[i][idir1]) for i in range(len(mass))])/len(mass),sum([(pos[i][idir2]) for i in range(len(mass))])/len(mass)]  # unweighted centroid 
    dpos1 = np.asarray([pos[i][idir1]-cent[0] for i in range(len(mass))])
    dpos2 = np.asarray([pos[i][idir2]-cent[1] for i in range(len(mass))])
    r2 = dpos1*dpos1+dpos2*dpos2  # radii of SF-gas from centroid 
    sortindex = np.argsort(r2)
    r2sort = r2[sortindex]
    msum = np.cumsum(mass[sortindex])
    for i in range(len(mass)):
        if msum[i] > frac*mtot:
            if r2sort[i] > 100*100: rh = mygals[igal].radii['stellar_half_mass']
            else: rh = np.sqrt(r2sort[i])
            break
    return rh,cent

def tage_init(cosmo,redshift):
    thubble = cosmo.age(redshift).value
    lz1 = np.arange(np.log10(1+redshift),np.log10(51),0.0001)
    tlz1 = thubble-cosmo.age(10**lz1-1.).value
    return lz1,tlz1

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

    # compute rhalf
        rhalf = np.zeros([3,len(mygals)])
        rhalfmass = np.zeros([3,len(mygals)])
        if conc:
                r20 = np.zeros([3,len(mygals)])
                r80 = np.zeros([3,len(mygals)])
        if plotvar=='mstar':
            for igal in range(len(mygals)):
                mass = np.array([sm[k] for k in mygals[igal].slist])
                pos = np.array([sp[k] for k in mygals[igal].slist])
                for idir0 in range(3): rhalf[idir0][igal],cent = compute_rfrac(idir0,mass,pos, frac)
    	    if igal%(len(mygals)/20)==0: print('%d logM*= %.3f %.3f c= [%.3f,%.3f] rh2d= %.3f %.3f %.3f rh3d= %.3f'%(igal,np.log10(ms[igal]),np.log10(sum(mass)),cent[0],cent[1],rhalf[0][igal],rhalf[1][igal],rhalf[2][igal],mygals[igal].radii['stellar_half_mass']))        
            with h5py.File(save_file, 'w') as f:
                hf.create_dataset('halfmass', data=np.array(rhalf))

        elif plotvar=='sfr':
            for igal in range(len(mygals)):
                mass = np.array([gsfr[k] for k in mygals[igal].glist])  # "mass" (for weighting) is actually now gas SFR
                pos = np.array([gp[k] for k in mygals[igal].glist])
                if sum(mass)==0: 
                    rhalf[0][igal] = rhalf[1][igal] = rhalf[2][igal] = 0
                    continue
                for idir0 in range(3): rhalf[idir0][igal],cent = compute_rfrac(idir0,mass,pos,frac)
    	    if igal%(len(mygals)/20)==0: print('%d logM*= %.3f c= [%.3f,%.3f] rh2d= %.3f %.3f %.3f rh3d= %.3f'%(igal,np.log10(ms[igal]),cent[0],cent[1],rhalf[0][igal],rhalf[1][igal],rhalf[2][igal],mygals[igal].radii['stellar_half_mass']))

        else:
            pgsnap = pg.Snap(snapfile,physical=True)
            print 'Computing',plotvar,'magnitudes for',len(pgsnap.stars),'stars.'
            allL = pg.snapshot.get_luminosities(pgsnap.stars,band=plotvar)
            print allL
            for igal in range(len(mygals)):
                mass = np.array([sm[k] for k in mygals[igal].slist])
                lum = np.array([allL[k] for k in mygals[igal].slist])
                #if igal%(len(mygals)/20)==0: print igal,mags[:3],lum[:3],min(mags),max(mags),min(lum),max(lum)
                pos = np.array([sp[k] for k in mygals[igal].slist])
                if not conc:
                    for idir0 in range(3): rhalf[idir0][igal],cent = compute_rfrac(idir0,lum,pos)
                    for idir0 in range(3): rhalfmass[idir0][igal],cent = compute_rfrac(idir0,mass,pos)
                    if igal%(len(mygals)/20)==0: print('%d logM*= %.3f c= [%.3f,%.3f] rh2d= %.2f %.2f %.2f rh2dMs= %.2f %.2f %.2f rh3d= %.3f'%(igal,np.log10(ms[igal]),cent[0],cent[1],rhalf[0][igal],rhalf[1][igal],rhalf[2][igal],rhalfmass[0][igal],rhalfmass[1][igal],rhalfmass[2][igal],mygals[igal].radii['stellar_half_mass']))
                elif conc:
                    for idir0 in range(3): r20[idir0][igal],cent = compute_rfrac(idir0,lum,pos, frac=0.2)
                    for idir0 in range(3): r80[idir0][igal],cent = compute_rfrac(idir0,lum,pos, frac=0.8)
                    if igal%(len(mygals)/20)==0: print('%d logM*= %.3f c= [%.3f,%.3f] r20= %.2f %.2f %.2f r80= %.2f %.2f %.2f rh3d= %.3f'%(igal,np.log10(ms[igal]),cent[0],cent[1],r20[0][igal],r20[1][igal],r20[2][igal],r80[0][igal],r80[1][igal],r80[2][igal],mygals[igal].radii['stellar_half_mass']))
            
            with h5py.File(save_file, 'a') as hf:
                if not conc:
                    hf.create_dataset(MODEL+'_'+WIND[iwind]+'_'+str(SNAP)+'_halflight', data=np.array(rhalf))
                    hf.create_dataset(MODEL+'_'+WIND[iwind]+'_'+str(SNAP)+'_halfmass', data=np.array(rhalfmass))
                elif conc:
                    hf.create_dataset(MODEL+'_'+WIND[iwind]+'_'+str(SNAP)+'_r20', data=np.array(r20))
                    hf.create_dataset(MODEL+'_'+WIND[iwind]+'_'+str(SNAP)+'_r80', data=np.array(r80))
