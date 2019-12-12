
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

def absmag_to_lum(absmags):
    # assuming absolute AB magnitudes 
    magab0 = 3631e-26 # W m^-2 Hz^-1
    pc = 3.085e16 # m
    
    flux = np.ones_like(absmags)*10.
    flux  = np.power(flux, absmags*-0.4)
    
    return flux*magab0*4.*np.pi*(10*pc)**2.


def compute_rfrac(idir0,mass,pos, frac=0.5):
    idir1 = (idir0+1)%3
    idir2 = (idir0+2)%3
    mtot = sum(mass)
    
    cent = [sum([(pos[i][idir1]) for i in range(len(mass))])/len(mass),sum([(pos[i][idir2]) for i in range(len(mass))])/len(mass)]  # unweighted centroid 
    dpos1 = np.asarray([pos[i][idir1]-cent[0] for i in range(len(mass))])
    dpos2 = np.asarray([pos[i][idir2]-cent[1] for i in range(len(mass))])
    
    r2 = dpos1*dpos1+dpos2*dpos2  # radii of SF-gas from centroid 
    sortindex = np.argsort(r2)
    r2sort = r2[sortindex]
    
    msum = np.cumsum(mass[sortindex])

    for i in range(len(mass)):
        if msum[i] > frac*mtot:
            if r2sort[i] > 100*100: 
                rh = mygals[igal].radii['stellar_half_mass']
            else: rh = np.sqrt(r2sort[i])
            break
    return rh,cent

# compile magnitudes in desired bands from pyloser file. 
# mtype choices are absmag/appmag/absmag_nodust/appmag_nodust.  also can be spec, iobjs, A_V, or L_FIR (in units of Lsun), in which case bandlist is not needed.
# bandlist is a list of the names of the desired bands, e.g. ['sdss_r','irac_1,'v'].  these must match the band names in the pyloser file.
def get_mags(loserfile,mtype,bandlist=None,verbose=False):
    hf = h5py.File(loserfile,'r')
    # read in pyloser parameters
    p = [i for i in hf.attrs.items()]
    params = dict([(p[i][0],p[i][1]) for i in range(len(p))])
    if verbose:
        print('Loaded params into dict:',params)
        print('hdf5 file keys: %s' % hf.keys())
    # get desired quantities
    for i in hf.keys():
        if i==mtype: mydata = list(hf[i])
    hf.close()

    # if A_V or L_FIR is requested, these can be directly returned
    if mtype == 'A_V' or mtype == 'L_FIR' or mtype == 'spec' or mtype == 'iobjs': return params,mydata

    # if specific bands are requested, find them in the pyloser file
    bands = params['bands']
    iband1 = iband2 = -1
    magindex = []
    for j in range(len(bandlist)):
        for i in range(len(bands)):
            if bands[i] == bandlist[j]: magindex.append(i)
            #if bands[i].decode('utf-8') == bandlist[j]: magindex.append(i)
    if len(magindex) != len(bandlist):
        print('Bands [',bandlist,'] not all found in pyloser file; available bands=',bands)
        exit()
    mags = []
    # compile desired band magnitudes
    for j in range(len(bandlist)):
        mags.append(np.asarray([m[magindex[j]] for m in mydata]))
    return params,mags

if __name__ == '__main__':
    
    save_file = '/home/sapple/simba_sizes/sizes/data/pyloser_sdss_r.h5'
    mtype = 'absmag'
    band = 'sdss_r'

    # define input file
    MODEL = 'm50n512'
    WIND = ['s50j7k', 's50noagn', 's50nojet', 's50nox']
    WIND = ['s50']
    mlim = 8.7
    mmax = 12.7
    snaps = ['151']

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
            params, mags = get_mags(loserfile, mtype, band, verbose=True)
            allL = absmag_to_lum(mags[0])
            
            for igal in range(len(mygals)):
                lum = np.array([allL[k] for k in mygals[igal].slist])
                pos = np.array([sp[k] for k in mygals[igal].slist])
                for idir0 in range(3): 
                    rhalf[idir0][igal],cent = compute_rfrac(idir0,lum,pos)
                    if igal%(len(mygals)/20)==0: print('%d logM*= %.3f c= [%.3f,%.3f] rh2d= %.2f %.2f %.2f rh3d= %.3f'%(igal,np.log10(ms[igal]),cent[0],cent[1],rhalf[0][igal],rhalf[1][igal],rhalf[2][igal],mygals[igal].radii['stellar_half_mass']))
                    
            with h5py.File(save_file, 'a') as hf:
                hf.create_dataset(MODEL+'_'+WIND[iwind]+'_'+str(SNAP)+'_halflight', data=np.array(rhalf))
                hf.create_dataset(MODEL+'_'+WIND[iwind]+'_'+str(SNAP)+'_halfmass', data=np.array(rhalfmass))
