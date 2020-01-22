import matplotlib.pyplot as plt
import h5py
import numpy as np
import os
import caesar
import sys
sys.path.append('/home/sapple/tools')
import plotmedian as pm

from size_plotting_methods import plot_data



rhalf_file = '/home/sapple/simba_sizes/sizes/data/halfradius_agn_R.h5'
plot_dir = '/home/sapple/simba_sizes/sizes/plots/'

winds = ['s50j7k', 's50nojet', 's50nox']
lines = ['-', '--', ':']
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
    ax.plot(bin_cent,ymean,lines[i],lw=2,color='c', label=wind)
    bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(np.log10(ms[ssfr<ssfrlim]),np.log10(rhalf[ssfr<ssfrlim]))
    ax.plot(bin_cent,ymean,lines[i],lw=2,color='m')

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
