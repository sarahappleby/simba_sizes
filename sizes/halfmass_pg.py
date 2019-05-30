
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


if len(sys.argv) < 3:
    print 'usage: MODEL WIND1 WIND2 ...'
    exit()

plotvar = 'mstar'
plotvar = 'R'  # options are 'bolometric', 'U', 'B', 'V', 'R', and 'K'
if len(sys.argv)==5: plotvar=sys.argv[4]

# define input file
MODEL = sys.argv[1]
#SNAP = int(sys.argv[2])
WIND = [sys.argv[2]]
#plt.rc('text', usetex=True)
mlim = 8.7
mmax = 12.7

snaps = ['151', '145', '125', '105', '090', '078', '062']

save_file = './halfradius.h5'
#plot_dir = '/home/sapple/simba_sizes/sizes/plots/'
#address = "21/3 viewforth gardens"

'''
if plotvar not in ['mstar','sfr']:
    try: sdss_bands = fsps.find_filter(plotvar)
    except: sys.exit('Filter %s not found'%plotvar)
    print('Doing rhalf in band %s, generating FSPS stellar pop...'%plotvar)
    spop = fsps.StellarPopulation(zcontinuous=1, sfh=0, logzsol=0.0, dust_type=2, dust2=0.2)
    #print spop.ssp_ages
'''

def plot_data(z):
  if z < 0.5: # Zhang+17 SDSS
    ms_data = np.linspace(9,12.0,20)
    alpha = 0.24
    beta = 1.33
    gamma = 10.17
    M0 = 6.49e11/0.68**2
    logRe_spiral = np.log10(gamma * (10**ms_data/M0)**alpha * (1+10**ms_data/M0)**(beta-alpha))
    plt.plot(ms_data,logRe_spiral,':',color='b',lw=3,label='SDSS-Blue')
    alpha = 0.17
    beta = 0.58
    gamma = 2.24
    M0 = 2.11e10/0.68**2
    logRe_etg = np.log10(gamma * (10**ms_data/M0)**alpha * (1+10**ms_data/M0)**(beta-alpha))
    plt.plot(ms_data,logRe_etg,':',color='r',lw=3,label='SDSS-Red')
    # Kelvin+11 GAMA
    #logRe_spiral = 0.4+0.3*(ms_data-9)
    #plt.plot(ms_data,logRe_spiral,'-.',color='b',lw=2,label='GAMA-Blue')
  if z >= 0.5 and z < 1.5: # vdWel+14 CANDELS
    logms = [9.25,9.75,10.25,10.75,11.25]
    logRe = [0.37,0.48,0.57,0.67,0.82]
    eRelo = [0.26,0.27,0.24,0.20,0.20]
    eRehi = [0.23,0.21,0.20,0.24,0.14]
    #plt.plot(logms,logRe,'o',color='k',ms=6,label='CANDELS LTG (van der Wel+14)')
    plt.errorbar(np.array(logms)-0.05,logRe,lw=1,yerr=[eRelo,eRehi],fmt='bo',ecolor='b',label='CANDELS LTG')
    logms = [9.25,9.75,10.25,10.75,11.25]
    logRe = [0.23,0.20,0.16,0.38,0.70]
    eRelo = [0.25,0.33,0.26,0.23,0.27]
    eRehi = [0.20,0.34,0.27,0.24,0.23]
    #plt.plot(logms,logRe,'o',color='k',ms=6,label='CANDELS ETG (van der Wel+14)')
    plt.errorbar(np.array(logms)+0.05,logRe,lw=1,yerr=[eRelo,eRehi],fmt='ro',ecolor='r',label='CANDELS ETG')
  if z >= 1.5 and z < 2.5:
    # Alcorn+16 CANDELS+ZFOURGE
    ms_data = np.linspace(9,11.5,5)
    re_data = 0.2*(ms_data-10)+0.4
    plt.plot(ms_data,re_data,'-',color='k',lw=4,label='CANDELS+ZFOURGE (Allen+16)')
    # vdWel+14 CANDELS
    logms = [9.75,10.25,10.75,11.25]
    logRe = [0.39,0.44,0.47,0.62]
    eRelo = [0.27,0.32,0.39,0.30]
    eRehi = [0.23,0.21,0.24,0.21]
    plt.plot(logms,logRe,'o',color='k',ms=6,label='CANDELS (van der Wel+14)')
    plt.errorbar(logms,logRe,lw=2,yerr=[eRelo,eRehi],color='k')

def compute_rhalf(idir0,mass,pos):
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
        if msum[i] > 0.5*mtot:
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
        if plotvar=='mstar':
            for igal in range(len(mygals)):
                mass = np.array([sm[k] for k in mygals[igal].slist])
                pos = np.array([sp[k] for k in mygals[igal].slist])
                for idir0 in range(3): rhalf[idir0][igal],cent = compute_rhalf(idir0,mass,pos)
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
                for idir0 in range(3): rhalf[idir0][igal],cent = compute_rhalf(idir0,mass,pos)
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
                for idir0 in range(3): rhalf[idir0][igal],cent = compute_rhalf(idir0,lum,pos)
                for idir0 in range(3): rhalfmass[idir0][igal],cent = compute_rhalf(idir0,mass,pos)
    	    if igal%(len(mygals)/20)==0: print('%d logM*= %.3f c= [%.3f,%.3f] rh2d= %.2f %.2f %.2f rh2dMs= %.2f %.2f %.2f rh3d= %.3f'%(igal,np.log10(ms[igal]),cent[0],cent[1],rhalf[0][igal],rhalf[1][igal],rhalf[2][igal],rhalfmass[0][igal],rhalfmass[1][igal],rhalfmass[2][igal],mygals[igal].radii['stellar_half_mass']))
            with h5py.File(save_file, 'a') as hf:
                hf.create_dataset(MODEL+'_'+WIND[0]+'_'+str(SNAP)+'_halflight', data=np.array(rhalf))
                hf.create_dataset(MODEL+'_'+WIND[0]+'_'+str(SNAP)+'_halfmass', data=np.array(rhalfmass))


"""
    #print len(ms),len(sigv3d),len(sigmastar)
    ssfr = 1.e9*sfr/ms
    ssfr = np.log10(ssfr+10**(-2.9+0.3*redshift))
    ssfrlim = -1.8+0.3*redshift
    rad = (rhalf[0]+rhalf[1]+rhalf[2])/3

    condition = (central)
    massbin,cvecbin,ebinlo,ebinhi = pm.runningmedian(np.log10(ms[condition]),ssfr[condition])
    cvec = ssfr# - np.interp(np.log10(ms),massbin,cvecbin)

    xvec = np.log10(ms[rad>0])
    yvec = np.log10(rad[rad>0])
    cvec = cvec[rad>0]
    pixsize = 1*(xvec-min(xvec))+0.5
    im = ax.scatter(xvec, yvec, c=cvec, s=pixsize, lw=0, cmap=plt.cm.jet_r)
    #plt.plot(xvec, yvec, 'o', c='k', ms=2)
    fig.colorbar(im,ax=ax,label=r'$\log$ sSFR (Gyr$^{-1}$)')
    if redshift <= 1.5:
        bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(np.log10(ms[ssfr>ssfrlim]),np.log10(rad[ssfr>ssfrlim]),bins=16)
        ax.plot(bin_cent,ymean,'-',lw=3,color='c',label='Simba-SF')
        #ax.errorbar(bin_cent,ymean,yerr=[ysiglo,ysighi],fmt='ro',c='c')
        bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(np.log10(ms[ssfr<ssfrlim]),np.log10(rad[ssfr<ssfrlim]),bins=16)
        ax.plot(bin_cent,ymean,'-',lw=3,color='m',label='Simba-Q')
        #ax.errorbar(bin_cent,ymean,yerr=[ysiglo,ysighi],fmt='ro',c='r')
    else:
        bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(xvec,yvec)
        ax.plot(bin_cent,ymean,'--',lw=3,color='g')
    #im.set_clim(-0.5,0.5)

plot_data(redshift)
#plot_data(0,'.')

plt.annotate('z=%g'%(np.round(redshift,1)), xy=(0.1, 0.9), xycoords='axes fraction',size=16,bbox=dict(boxstyle="round", fc="w"))

plt.minorticks_on()
plt.xlim(mlim+0.3,mmax)
plt.ylim(-0.2,1.8-0.3*redshift)
plt.xlabel(r'$\log\ M_{*}$',fontsize=20)
plt.ylabel(r'$\log\ R_{half}$' ,fontsize=20)
plt.legend(loc='lower right')

plt.savefig(plot_dir+'halfmass_%s.pdf'%MODEL, bbox_inches='tight', format='pdf')

plt.show()
"""
