import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import os
import h5py
import caesar
import numpy as np

import sys
sys.path.append('/home/sapple/tools')
import plotmedian as pm

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
                                                        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def plot_data(ax, z):
        if z < 0.5: # Zhang+17 SDSS
                ms_data = np.linspace(9,12.0,20)
                alpha = 0.24
                beta = 1.33
                gamma = 10.17
                M0 = 6.49e11/0.68**2
                logRe_spiral = np.log10(gamma * (10**ms_data/M0)**alpha * (1+10**ms_data/M0)**(beta-alpha))
                ax.plot(ms_data,logRe_spiral,':',color='b',lw=6,label='SDSS-LTG (Zhang+17)', markersize=3)
                alpha = 0.17
                beta = 0.58
                gamma = 2.24
                M0 = 2.11e10/0.68**2
                logRe_etg = np.log10(gamma * (10**ms_data/M0)**alpha * (1+10**ms_data/M0)**(beta-alpha))
                ax.plot(ms_data,logRe_etg,':',color='r',lw=6,label='SDSS-ETG (Zhang+17)', markersize=3)
                # Kelvin+11 GAMA
                #logRe_spiral = 0.4+0.3*(ms_data-9)
                #plt.plot(ms_data,logRe_spiral,'-',color='k',lw=4,label='GAMA (Baldry+12)')
        if z >= 0.5 and z < 1.5: # vdWel+14 CANDELS
                logms = [9.25,9.75,10.25,10.75,11.25]
                logRe = [0.37,0.48,0.57,0.67,0.82]
                eRelo = [0.26,0.27,0.24,0.20,0.20]
                eRehi = [0.23,0.21,0.20,0.24,0.14]
                #plt.plot(logms,logRe,'o',color='k',ms=6,label='CANDELS LTG (van der Wel+14)')
                ax.errorbar(np.array(logms)-0.05,logRe,lw=3,yerr=[eRelo,eRehi],fmt='bo',ecolor='b',label='CANDELS LTG', markersize=6)
                logms = [9.25,9.75,10.25,10.75,11.25]
                logRe = [0.23,0.20,0.16,0.38,0.70]
                eRelo = [0.25,0.33,0.26,0.23,0.27]
                eRehi = [0.20,0.34,0.27,0.24,0.23]
                #plt.plot(logms,logRe,'o',color='k',ms=6,label='CANDELS ETG (van der Wel+14)')
                ax.errorbar(np.array(logms)+0.05,logRe,lw=3,yerr=[eRelo,eRehi],fmt='ro',ecolor='r',label='CANDELS ETG', markersize=6)
        if z >= 1.5 and z < 2.5:
                # Alcorn+16 CANDELS+ZFOURGE
                ms_data = np.linspace(9,11.5,5)
                re_data = 0.2*(ms_data-10)+0.4
                ax.plot(ms_data,re_data,'-',color='k',lw=3,label='CANDELS+ZFOURGE (Allen+16)')
                # vdWel+14 CANDELS
                logms = [9.75,10.25,10.75,11.25]
                logRe = [0.39,0.44,0.47,0.62]
                eRelo = [0.27,0.32,0.39,0.30]
                eRehi = [0.23,0.21,0.24,0.21]
                ax.plot(logms,logRe,'o',color='k',label='CANDELS (van der Wel+14)', markersize=6)
                ax.errorbar(logms,logRe,lw=3,yerr=[eRelo,eRehi],color='k')


rhalf_file = '/home/sapple/simba_sizes/sizes/data/halfradius.h5'
plot_dir = '/home/sapple/simba_sizes/sizes/plots/'

simba_snaps = ['062', '078', '090', '105', '125', '151']
simba_z = [3.0, 2.0, 1.5, 1.0, 0.5, 0.0]

model = 'm100n1024'
wind = 's50j7k'

mmin = 8.7
mmax = 12.5

# one subplot per redshift

#fig, axs = plt.subplots(3, 2)
#axes = axs.reshape(-1)

cmap = plt.get_cmap('jet_r')
cmap = truncate_colormap(cmap, 0.05, 1.0)

fig = plt.figure(figsize=(15, 17))

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
                central = f['central'].value
                ms = f['stellar_mass'].value
                sfr = f['sfr'].value

	with h5py.File(rhalf_file, 'r') as r:
		rhalf = r[model+'_'+wind+'_'+simba_snaps[i]+'_halflight'].value
		#rhalfmass = r[model+'_'+wind+'_'+simba_snaps[i]+'_halfmass'].value
        rhalf = np.sum(rhalf, axis=0) / 3

	ssfr = 1.e9*sfr/ms
        #ssfr = np.log10(ssfr+10**(-2.9+0.3*simba_z[i]))
	ssfr = np.log10(ssfr)
	ssfrlim = -1.8+0.3*simba_z[i]

	# plot simba galaxy points
	massbin,cvecbin,ebinlo,ebinhi = pm.runningmedian(np.log10(ms[central]),ssfr[central])
        simba_x = np.log10(ms[rhalf>0])  
        simba_y = np.log10(rhalf[rhalf>0])
        simba_c = ssfr - np.interp(np.log10(ms),massbin,cvecbin)
        simba_c = ssfr
        simba_c[simba_c < -2.5] = -2.5

	pixsize = 1*(simba_x-min(simba_x))+0.5
	
	ax = fig.add_subplot(3, 2, i+1)

	im = ax.scatter(simba_x, simba_y, c=simba_c, s=pixsize, lw=0, cmap=cmap, label='Simba')
	cbar = fig.colorbar(im,ax=ax, label=r'$\textrm{log} (\textrm{sSFR} / \textrm{Gyr}^{-1})$')
	cbar.ax.tick_params(labelsize=10)

	# plot simba red and blue/ all galaxies
	if simba_z[i] <= 1.5:
        	bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(np.log10(ms[ssfr>ssfrlim]),np.log10(rhalf[ssfr>ssfrlim]))
                ax.plot(bin_cent,ymean,'--',lw=4,color='c', label='Star forming')   
		bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(np.log10(ms[ssfr<ssfrlim]),np.log10(rhalf[ssfr<ssfrlim]))
                ax.plot(bin_cent,ymean,'--',lw=4,color='m', label='Passive')
	else:
		bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(simba_x,simba_y)
		ax.plot(bin_cent,ymean,'--',lw=4,color='b', label='Simba')
	
	# plot redshift observational data
	plot_data(ax, simba_z[i])

	ax.annotate('z=%g'%(np.round(simba_z[i],1)), xy=(0.1, 0.9), xycoords='axes fraction',size=16,bbox=dict(boxstyle="round", fc="w"))
	ax.minorticks_on()
	ax.set_xlim(mmin,mmax)
	ax.set_ylim(-0.3,1.8-0.3*simba_z[i])
	ax.set_xlabel(r'$\log\ (M_{*} / M_{\odot})$',fontsize=16)
	ax.set_ylabel(r'$\log\ (R_{half} / \textrm{kpc})$' ,fontsize=16)
	ax.legend(loc='lower right', fontsize=8)
	
plt.savefig(plot_dir+'halflight_'+model+'.png', bbox_inches='tight', format='png')
plt.clf()

