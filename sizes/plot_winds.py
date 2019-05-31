import matplotlib.pyplot as plt
import h5py
import numpy as np
import os
import caesar
import sys
sys.path.append('/home/sapple/tools')
import plotmedian as pm


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


rhalf_file = '/home/sapple/simba_sizes/sizes/data/halfradius_agn.h5'
plot_dir = '/home/sapple/simba_sizes/sizes/plots/'

winds = ['s50j7k', 's50nojet', 's50nox']
lines = ['--', '-.', ':']
snap = '151'
model = 'm50n512'
simba_z = 0.

mmin = 8.7
mmax = 12.5

fig, ax = plt.subplots()
for i, wind in enumerate(winds):

    caesar_data = '/home/sapple/simba_sizes/sizes/data/'+model+'_'+wind+'_'+snap+'_caesar_data.h5'
    if not os.path.isfile(caesar_data):
        infile = '/home/rad/data/'+model+'/'+wind+'/Groups/'+model+'_'+snap+'.hdf5'
        sim = caesar.load(infile, LoadHalo=False)
        central = np.asarray([j.central for j in sim.galaxies])
        ms = np.asarray([j.masses['stellar'].in_units('Msun') for j in sim.galaxies])
        sfr = np.array([j.sfr.in_units('Msun/yr') for j in sim.galaxies])
        sfr[np.where(sfr == 1.)[0]] = 0.
        
        with h5py.File(caesar_data, 'w') as f:
            f.create_dataset('central', data=np.array(central))
            f.create_dataset('stellar_mass', data=np.array(ms))
            f.create_dataset('sfr', data=np.array(sfr))
    else:
        with h5py.File(caesar_data, 'r') as f:
            central = f['central'].value
            ms = f['stellar_mass'].value
            sfr = f['sfr'].value

    ssfr = np.log10(1.e9*sfr/ms)
    ssfrlim = -1.8+0.3*simba_z

    with h5py.File(rhalf_file, 'r') as r:
        rhalf = r[model+'_'+wind+'_'+snap+'_halflight'].value
    rhalf = np.sum(rhalf, axis=0) / 3.
    
    bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(np.log10(ms[ssfr>ssfrlim]),np.log10(rhalf[ssfr>ssfrlim]))
    ax.plot(bin_cent,ymean,lines[i],lw=3,color='c', label=wind)
    bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(np.log10(ms[ssfr<ssfrlim]),np.log10(rhalf[ssfr<ssfrlim]))
    ax.plot(bin_cent,ymean,lines[i],lw=3,color='m')

plot_data(simba_z)
plt.annotate('z=%g'%(np.round(simba_z,1)), xy=(0.1, 0.9), xycoords='axes fraction',size=16,bbox=dict(boxstyle="round", fc="w"))
plt.minorticks_on()
#plt.xlim(mmin+0.3,mmax)
#plt.ylim(-0.2,1.8-0.3*simba_z)
plt.xlabel(r'$\log\ M_{*}$',fontsize=20)
plt.ylabel(r'$\log\ R_{half}$' ,fontsize=20)
plt.legend(loc='lower right')

plt.savefig(plot_dir+'halfmass_agn.png', bbox_inches='tight', format='png')
plt.clf()
