# SA : include method to read in and plot data for m25n512



import matplotlib.pyplot as plt
import matplotlib.colors as colors
from mpl_toolkits.axes_grid1 import make_axes_locatable

import os
import h5py
import caesar
import numpy as np

import sys
sys.path.append('/home/sapple/tools')
import plotmedian as pm

from size_plotting_methods import *

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=14)

if __name__ == '__main__':

	plot_dir = '/home/sapple/simba_sizes/sizes/plots/'

	simba_snaps = ['078', '105', '151']
	simba_z = [2.0, 1.0, 0.0]

	mtype = 'abs'

	model = 'm100n1024'
	wind = 's50'

	mmin = 10.
	mmax = 12.5

	cmap = plt.get_cmap('jet_r')
	cmap = truncate_colormap(cmap, 0.05, 1.0, alpha=0.5)

	fig, ax = plt.subplots(1, 3, sharey=True, figsize=(15, 6))

	for i in range(len(simba_z)):

		caesar_data = '/home/sapple/simba_sizes/sizes/data/'+model+'_'+wind+'_'+simba_snaps[i]+'_caesar_data.h5'
		if not os.path.isfile(caesar_data):
			infile = '/home/rad/data/'+model+'/'+wind+'/Groups/'+model+'_'+simba_snaps[i]+'.hdf5'
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

		if simba_z[i] == 0.:
			rhalf_file = '/home/sapple/simba_sizes/sizes/data/pyloser_sdss_r.h5'
		else:
			rhalf_file = '/home/sapple/simba_sizes/sizes/data/pyloser_v.h5'
		with h5py.File(rhalf_file, 'r') as r:
			rhalf = r[mtype+'_'+model+'_'+wind+'_'+simba_snaps[i]][:]
		rhalf = np.sum(rhalf, axis=0) / 3
		ssfr = 1.e9*sfr/ms
		ssfr = np.log10(ssfr)
		ssfrlim = -1.8+0.3*simba_z[i]

		# plot simba galaxy points
		massbin,cvecbin,ebinlo,ebinhi = pm.runningmedian(np.log10(ms[central]),ssfr[central])
		simba_x = np.log10(ms[(rhalf>0) & (ms > 10**mmin)])  
		simba_y = np.log10(rhalf[(rhalf>0) & (ms > 10**mmin)])
		simba_c = ssfr[(rhalf>0) & (ms > 10**mmin)]
		simba_c[simba_c < -2.5] = -2.5
		simba_c[simba_c > 1.] = 1.

		pixsize = 1*(simba_x-min(simba_x))+0.5

		im = ax[i].scatter(simba_x, simba_y, c=simba_c, s=2, lw=0, cmap=cmap)
		#im = ax[i].hexbin(simba_x, simba_y, C=simba_c, cmap=cmap, gridsize=30)

		# plot redshift observational data
		plot_data(ax[i], simba_z[i])

		# plot high resolution (m25n512) data
		plot_high_res(ax[i], simba_snaps[i], mtype)

		# plot simba red and blue/ all galaxies
		if simba_z[i] <= 2.5:
			bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(np.log10(ms[ssfr>ssfrlim]),np.log10(rhalf[ssfr>ssfrlim]))
			ax[i].plot(bin_cent,ymean,'-',lw=3.,color='dodgerblue', label='Star forming')   
			bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(np.log10(ms[ssfr<ssfrlim]),np.log10(rhalf[ssfr<ssfrlim]))
			ax[i].plot(bin_cent,ymean,'-',lw=3.,color='m', label='Passive')
		else:
			bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(simba_x,simba_y)
			ax[i].plot(bin_cent,ymean,'-.',marker='o', lw=3.,color='dodgerblue', label='Simba')


		res_line = np.log10(2.8*0.5 / ((1+simba_z[i])*0.68))
		ax[i].annotate('z=%g'%(np.round(simba_z[i],1)), xy=(0.05, 0.93), xycoords='axes fraction',size=16,)
		ax[i].minorticks_on()
		ax[i].set_xlim(mmin,mmax)
		ax[i].axhline(res_line, linestyle='--', color='k', lw=1.5, )
		ax[i].set_ylim(-0.2, 1.4)
		ax[i].set_xlabel(r'$\log\ (M_{*} / M_{\odot})$',fontsize=15)
		if i == 0:
			loc = 1
		else:
			loc = 4
		do_legend_order(ax[i], simba_z[i], loc)

		if i ==0:
			ax[i].set_ylabel(r'$\log\ (R_{half} / \textrm{kpc})$' ,fontsize=15)
		elif i == 2:
			cbaxes = fig.add_axes([0.9, 0.11, 0.02, 0.77])
			cb = plt.colorbar(im, cax=cbaxes, label=r'$\textrm{log} (\textrm{sSFR} / \textrm{Gyr}^{-1})$', drawedges=False)
			#cbar = fig.colorbar(im,ax=ax.ravel().tolist(), label=r'$\textrm{log} (\textrm{sSFR} / \textrm{Gyr}^{-1})$')
			#cbar.ax.tick_params(labelsize=12)

	fig.subplots_adjust(left=0.08, right=0.88, wspace=0.1)
	plt.savefig(plot_dir+mtype+'_halflight_evol_'+model+'.png', bbox_inches='tight', format='png', pad_inches=0)
	plt.show()
