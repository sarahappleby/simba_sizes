import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import h5py
import numpy as np
import os
import caesar
import sys
sys.path.append('/home/sapple/tools')
import plotmedian as pm

from size_plotting_methods import plot_data, plot_sdss_sf

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=16)


rhalf_file = '/home/sapple/simba_sizes/sizes/data/pyloser_sdss_r.h5'
plot_dir = '/home/sapple/simba_sizes/sizes/plots/'

winds = ['s50', 's50nox', 's50nojet', 's50noagn']
lines = ['-', '--', '-.', ':']
snap = '151'
model = 'm50n512'
simba_z = 0.
mtype = 'abs'

mmin = 10.
mmax = 12.5


fig, ax = plt.subplots(figsize=(7, 7))

plot_data(ax, simba_z)

for i, w in enumerate(winds):
    ax.plot([0,1],[0,1],linestyle=lines[i], color='k', label=w)


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
            central = f['central'][:]
            ms = f['stellar_mass'][:]
            sfr = f['sfr'][:]

    ssfr = np.log10(1.e9*sfr/ms)
    ssfrlim = -1.8+0.3*simba_z

    with h5py.File(rhalf_file, 'r') as r:
        rhalf = r[mtype+'_'+model+'_'+wind+'_'+snap][:]
    rhalf = np.sum(rhalf, axis=0) / 3.
    
    bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(np.log10(ms[ssfr>ssfrlim]),np.log10(rhalf[ssfr>ssfrlim]))
    ax.plot(bin_cent,ymean,lines[i],lw=2,color='blue')
    bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(np.log10(ms[ssfr<ssfrlim]),np.log10(rhalf[ssfr<ssfrlim]))
    ax.plot(bin_cent,ymean,lines[i],lw=2,color='m')


res_line = np.log10(2.8*0.5 / ((1+simba_z)*0.68))
ax.axhline(res_line, linestyle='--', color='k', lw=1.5, )
plt.annotate('z=%g'%(np.round(simba_z,1)), xy=(0.05, 0.93), xycoords='axes fraction',size=15,bbox=dict(boxstyle="round", fc="w"))
plt.minorticks_on()
plt.xlim(mmin, mmax)
plt.ylim(0.2, )
plt.xlabel(r'$\log\ M_{*}$',fontsize=20)
plt.ylabel(r'$\log\ R_{half}$' ,fontsize=20)
plt.legend(loc='lower right', fontsize=14)

plt.savefig(plot_dir+'halflight_agn.png', bbox_inches='tight', format='png')
plt.show()
