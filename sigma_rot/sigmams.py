"""
From Romeel.
"""
import caesar
import sys
import pylab as plt
import numpy as np
from readgadget import *
import plotmedian as pm

if len(sys.argv) < 3:
    print 'usage: MODEL SNAP WIND1 WIND2 ...'
    exit()

# define input file
MODEL = sys.argv[1]
SNAP = int(sys.argv[2])
WIND = sys.argv[3:]
#plt.rc('text', usetex=True)
mlim = 8.7
mmax = 12.5

def plot_data(redshift):
  if redshift < 0.5:	
    infile = 'Observations/KormendyHo2013/KH13.dat'
    hubtype,mbh,mbhlo,mbhhi,sig,esig = np.loadtxt(infile,usecols=(2,11,12,13,14,15),unpack=True)
    #print mbh,mbhlo,mbhhi,sig,esig
    embh = [np.log10(mbh)-np.log10(mbh-mbhlo),np.log10(mbhhi+mbh)-np.log10(mbh)]
    elogsig = [np.log10(sig)-np.log10(sig-esig),np.log10(sig+esig)-np.log10(sig)]
    plt.errorbar(np.log10(sig),np.log10(mbh*1.e6),fmt='.',yerr=embh,xerr=elogsig,lw=1,color='grey',label='Kormendy+Ho')

# load in input file
fig,ax = plt.subplots()
for iwind in range(0,len(WIND)):
    snapfile = '/home/rad/data/%s/%s/snap_%s_%03d.hdf5' % (MODEL,WIND[iwind],MODEL,SNAP)
    infile = '/home/rad/data/%s/%s/Groups/%s_%03d.hdf5' % (MODEL,WIND[iwind],MODEL,SNAP)
    sim = caesar.load(infile)
    redshift = sim.simulation.redshift
    
    ids = np.asarray([i.GroupID for i in sim.galaxies if i.central == 1])
    ms = np.asarray([i.masses['stellar'] for i in sim.galaxies if i.central == 1])
    mbh = np.asarray([i.masses['bh'] for i in sim.galaxies if i.central == 1])
    sfr = np.asarray([i.sfr for i in sim.galaxies if i.central == 1])
    met = np.asarray([i.metallicities['sfr_weighted'] for i in sim.galaxies if i.central == 1])
    sigmastar = np.asarray([i.velocity_dispersions['stellar'] for i in sim.galaxies if i.central == 1])
    rad = np.asarray([i.radii['stellar_half_mass'] for i in sim.galaxies if i.central == 1])

    ms_sat = np.asarray([i.masses['stellar'] for i in sim.galaxies if i.central != 1])
    mbh_sat = np.asarray([i.masses['bh'] for i in sim.galaxies if i.central != 1])
    sigmastar_sat = np.asarray([i.velocity_dispersions['stellar'] for i in sim.galaxies if i.central != 1])
    rad_sat = np.asarray([i.radii['stellar_half_mass'] for i in sim.galaxies if i.central != 1])
  
    cents = np.asarray([i for i in sim.galaxies if i.central == 1])
    sv = readsnap(snapfile,'vel','star',units=1,suppress=1) # physical km/s
    sigv3d = []
    for g in cents:
        svgal = np.array([sv[k] for k in g.slist])
        svcent = np.mean(svgal,axis=0)
	svgal = (svgal-svcent)*(svgal-svcent)
	sigv = np.sqrt(sum(svgal)/len(svgal))
	sigv3d.append(np.sqrt(sigv[0]*sigv[0]+sigv[1]*sigv[1]+sigv[2]*sigv[2]))
	#print np.log10(g.masses['stellar']),sigv,sigv3d[:1],g.velocity_dispersions['stellar']

    #print len(ms),len(sigv3d),len(sigmastar)
    logms = np.log10(ms)
    logmbh = np.log10(mbh)
    logsig = np.log10(np.asarray(sigv3d))
    logms_sat = np.log10(ms_sat)
    logmbh_sat = np.log10(mbh_sat)
    logsig_sat = np.log10(sigmastar)
    Zmet = np.log10(met/0.0189+1.e-6)
    ssfr = 1.e9*sfr/ms
    ssfr = np.log10(ssfr+10**(-2.5+0.3*redshift))
    pixcolor = np.log10(rad)
    pixsize = 4*(np.log10(ms/min(ms))+1)

    cvec = np.log10(rad)
    cvec_sat = np.log10(rad_sat)
    massbin,cvecbin,ebinlo,ebinhi = pm.runningmedian(np.log10(ms[sfr>0]),np.log10(rad[sfr>0]))
    cvec = cvec - np.interp(np.log10(ms),massbin,cvecbin)
    cvec_sat = cvec_sat - np.interp(np.log10(ms_sat),massbin,cvecbin)

    if iwind==0: 
        #ax.plot(logms,logsig_sat, 'x', c='grey', ms=1, label='Sats')
        im = ax.scatter(logms,logsig, c=cvec, s=pixsize, lw=0, cmap=plt.cm.jet_r, label='Centrals')
        fig.colorbar(im,ax=ax,label=r'$\Delta r_*$ (kpc)')
        im.set_clim(-0.2,0.2)
    else: im = ax.scatter(logms_sat, logsig_sat, c='k', s=pixsize, lw=0)
    #plt.plot(logms,Zmet,ms=1,lw=0,c='r')
    #print logms,Zmet

#plot_data(0)
#plot_data(0,'.')

plt.annotate('z=%g,%s'%(np.round(redshift,1),WIND[0]), xy=(0.1, 0.9), xycoords='axes fraction',size=16,bbox=dict(boxstyle="round", fc="w"))

plt.minorticks_on()
plt.xlim(mlim,mmax)
plt.ylim(1.7,2.7)
plt.xlabel(r'$\log\ M_{*}$',fontsize=16)
plt.ylabel(r'$\log\ \sigma_{*}$' ,fontsize=16)
plt.legend(loc='lower right')

plt.savefig('sigmams.%s.pdf'%MODEL, bbox_inches='tight', format='pdf')

plt.show()

