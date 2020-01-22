import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

import os
import h5py
import caesar
import numpy as np

from size_plotting_methods import truncate_colormap, plot_data

import sys
sys.path.append('/home/sapple/tools')
import plotmedian as pm

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=16)

plot_dir = '/home/sapple/simba_sizes/sizes/plots/'

snap = '151'
mtype = 'abs'
model = 'm25n512'
wind = 's50'
simba_z = 0.

mmin = 10.
mmax = 12.5

cmap = plt.get_cmap('jet_r')
cmap = truncate_colormap(cmap, 0.03, 1.0)

fig, ax = plt.subplots(1, 1, sharey=True, figsize=(7, 7))

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

if simba_z == 0.:
    rhalf_file = '/home/sapple/simba_sizes/sizes/data/pyloser_sdss_r.h5'
else:
    rhalf_file = '/home/sapple/simba_sizes/sizes/data/pyloser_v.h5'
with h5py.File(rhalf_file, 'r') as r:
    rhalf = r[mtype+'_'+model+'_'+wind+'_'+snap].value
rhalf = np.sum(rhalf, axis=0) / 3
	
ssfr = 1.e9*sfr/ms
ssfr = np.log10(ssfr)
ssfrlim = -1.8+0.3*simba_z

# plot simba galaxy points
massbin,cvecbin,ebinlo,ebinhi = pm.runningmedian(np.log10(ms[central]),ssfr[central])
simba_x = np.log10(ms[rhalf>0])  
simba_y = np.log10(rhalf[rhalf>0])
simba_c = ssfr - np.interp(np.log10(ms),massbin,cvecbin)
simba_c = ssfr
simba_c[simba_c < -2.5] = -2.5
simba_c[simba_c > 1.] = 1.

pixsize = 1*(simba_x-min(simba_x))+0.5
	
im = ax.scatter(simba_x, simba_y, c=simba_c, s=pixsize, lw=0, cmap=cmap)
		
# plot simba red and blue/ all galaxies
if simba_z <= 2.5:
    bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(np.log10(ms[ssfr>ssfrlim]),np.log10(rhalf[ssfr>ssfrlim]))
    ax.plot(bin_cent,ymean,'--',lw=4,color='dodgerblue', label='Star forming')   
    bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(np.log10(ms[ssfr<ssfrlim]),np.log10(rhalf[ssfr<ssfrlim]))
    ax.plot(bin_cent,ymean,'--',lw=4,color='m', label='Passive')
else:
    bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(simba_x,simba_y)
    ax.plot(bin_cent,ymean,'--',lw=4,color='b', label='Simba')
	
# plot redshift observational data
plot_data(ax, simba_z)

res_line = np.log10(2.8*0.5 / ((1+simba_z)*0.68))
ax.annotate('z=%g'%(np.round(simba_z,1)), xy=(0.05, 0.9), xycoords='axes fraction',size=16,bbox=dict(boxstyle="round", fc="w"))
ax.minorticks_on()
#ax.set_xlim(mmin,mmax)
ax.axhline(res_line, linestyle='--', color='k', lw=1.5, )
ax.set_ylim(-0.3, 1.3)
ax.set_xlabel(r'$\log\ (M_{*} / M_{\odot})$',fontsize=16)
ax.legend(fontsize=11)

ax.set_ylabel(r'$\log\ (R_{half} / \textrm{kpc})$' ,fontsize=16)
cbar = fig.colorbar(im,ax=ax, label=r'$\textrm{log} (\textrm{sSFR} / \textrm{Gyr}^{-1})$')
cbar.ax.tick_params(labelsize=12)

plt.savefig(plot_dir+mtype+'_halflight_'+model+'_'+snap+'.png', bbox_inches='tight', format='png')
plt.show()
