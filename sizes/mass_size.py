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
        if z < 0.5: # Zhang+19 SDSS
                ms_data = np.linspace(9,12.0,20)
                alpha = 0.23
                beta = 0.41
                gamma = 10.72
                M0 = 4.9e11/0.68**2
                logRe_sf = np.log10(gamma * (10**ms_data/M0)**alpha * (1+10**ms_data/M0)**(beta-alpha))
                ax.plot(ms_data,logRe_sf,':',color='b',lw=6,label='SDSS-SF (Zhang+19)', markersize=3)
                alpha = -0.02
                beta = 0.65
                gamma = 1.58
                M0 = 1.11e10/0.68**2
                logRe_q = np.log10(gamma * (10**ms_data/M0)**alpha * (1+10**ms_data/M0)**(beta-alpha))
                ax.plot(ms_data,logRe_q,':',color='r',lw=6,label='SDSS-Q (Zhang+19)', markersize=3)
                # Kelvin+11 GAMA
                #logRe_spiral = 0.4+0.3*(ms_data-9)
                #plt.plot(ms_data,logRe_spiral,'-',color='k',lw=4,label='GAMA (Baldry+12)')
        if z >= 0.5 and z < 1.5: # vdWel+14 CANDELS
                #average of z = 0.75, 1.25, table 2
                logms = [9.25,9.75,10.25,10.75,11.25]
                logRe = [0.4, 0.52, 0.61, 0.71, 0.86]
                eRelo = [0.25, 0.25, 0.25, 0.22, 0.17]
                eRehi = [0.23, 0.21, 0.09, 0.16, 0.18]
                ax.errorbar(np.array(logms)-0.05,logRe,lw=3,yerr=[eRelo,eRehi],fmt='bo',ecolor='b',label='CANDELS LTG', markersize=6)
                logms = [9.25,9.75,10.25,10.75,11.25]
                logRe = [0.23, 0.195, 0.165, 0.375, 0.695]
                eRelo = [0.25, 0.33, 0.25, 0.21, 0.18]
                eRehi = [0.2, 0.24, 0.23, 0.22, 0.20]
                ax.errorbar(np.array(logms)+0.05,logRe,lw=3,yerr=[eRelo,eRehi],fmt='ro',ecolor='r',label='CANDELS ETG', markersize=6)
        if z >= 1.5 and z < 2.5:
                # Allen+17 CANDELS+ZFOURGE, table 1
                ms_data = np.linspace(9,11.5,5)
                re_data = 0.2*(ms_data-10)+0.4
                ax.plot(ms_data,re_data,'-',color='k',lw=3,label='CANDELS+ZFOURGE (Allen+17)')
                # vdWel+14 CANDELS, table 2
                # average of z = 1.75, 2.25, table 2
                logms = [9.75, 10.25, 10.75, 11.25]
                logRe = [0.39, 0.48, 0.57, 0.67]
                eRelo = [0.26, 0.26, 0.27, 0.21]
                eRehi = [0.22, 0.2, 0.18, 0.19]
                ax.errorbar(np.array(logms)-0.05,logRe,lw=3,yerr=[eRelo,eRehi],fmt='bo',ecolor='b',label='CANDELS LTG', markersize=6)
                logms = [9.75, 10.25, 10.75, 11.25]
                logRe = [0.22, -0.01, 0.14, 0.41]
                eRelo = [0.24, 0.31, 0.26, 0.19]
                eRehi = [0.26, 0.36, 0.38, 0.24]
                ax.errorbar(np.array(logms)+0.05,logRe,lw=3,yerr=[eRelo,eRehi],fmt='ro',ecolor='r',label='CANDELS ETG', markersize=6) 
                """ ltg + etg, table 5
                logms = [9.75,10.25,10.75,11.25]
                logRe = [0.39,0.44,0.47,0.62]
                eRelo = [0.27,0.32,0.39,0.30]
                eRehi = [0.23,0.21,0.24,0.21]
                ax.plot(logms,logRe,'o',color='k',label='CANDELS (van der Wel+14)', markersize=6)
                ax.errorbar(logms,logRe,lw=3,yerr=[eRelo,eRehi],color='k')
                """

rhalf_file = '/home/sapple/simba_sizes/sizes/data/halfradius.h5'
plot_dir = '/home/sapple/simba_sizes/sizes/plots/'

simba_snaps = ['078', '105', '151']
simba_z = [2.0, 1.0, 0.0]

model = 'm100n1024'
wind = 's50j7k'

mmin = 9.5
mmax = 12.5

# one subplot per redshift

#fig, axs = plt.subplots(3, 2)
#axes = axs.reshape(-1)

cmap = plt.get_cmap('jet_r')
cmap = truncate_colormap(cmap, 0.05, 1.0)

fig = plt.figure(figsize=(24, 6))

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
		#rhalf = r[model+'_'+wind+'_'+simba_snaps[i]+'_halfmass'].value
        rhalf = np.sum(rhalf, axis=0) / 3
        #rhalf = (rhalf[0]**2 + rhalf[1]**2 + rhalf[2]**2)**0.5
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
        simba_c[simba_c > 1.] = 1.

	pixsize = 1*(simba_x-min(simba_x))+0.5
	
	ax = fig.add_subplot(1, 3, i+1)

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

        res_line = np.log10(2.8*0.5 / ((1+simba_z[i])*0.68))
	ax.annotate('z=%g'%(np.round(simba_z[i],1)), xy=(0.1, 0.9), xycoords='axes fraction',size=16,bbox=dict(boxstyle="round", fc="w"))
	ax.minorticks_on()
	ax.set_xlim(mmin,mmax)
        ax.axhline(res_line, linestyle='--', color='k', lw=1.5, )
	#ax.set_ylim(-0.3,1.8-0.3*simba_z[i])
        ax.set_ylim(-0.25, 1.3)
	ax.set_xlabel(r'$\log\ (M_{*} / M_{\odot})$',fontsize=16)
	ax.set_ylabel(r'$\log\ (R_{half} / \textrm{kpc})$' ,fontsize=16)
	ax.legend(loc='lower right', fontsize=9)
	
plt.savefig(plot_dir+'halflight_'+model+'.png', bbox_inches='tight', format='png')
plt.clf()

