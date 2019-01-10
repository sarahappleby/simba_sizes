import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

import os
import h5py
import caesar
import numpy as np

import sys
sys.path.append('/home/sapple/tools')
import plotmedian as pm



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
                ax.plot(ms_data,re_data,'-',color='k',lw=3,label='CANDELS+ZFOURGE (Alcorn+16)')
                # vdWel+14 CANDELS
                logms = [9.75,10.25,10.75,11.25]
                logRe = [0.39,0.44,0.47,0.62]
                eRelo = [0.27,0.32,0.39,0.30]
                eRehi = [0.23,0.21,0.24,0.21]
                ax.plot(logms,logRe,'o',color='k',label='CANDELS (van der Wel+14)', markersize=6)
                ax.errorbar(logms,logRe,lw=3,yerr=[eRelo,eRehi],color='k')


data_dir = '/home/sapple/simba_sizes/data/'
plot_dir = '/home/sapple/simba_sizes/plots/'

mufasa_snaps = ['070', '085', '095', '105', '125', '126', '135']
simba_snaps = ['062', '078', '090', '105', '125', '145', '151']

mufasa_z = [3.0, 2.0, 1.5, 1.0, 0.5, 0.0]
simba_z = [3.0, 2.0, 1.5, 1.0, 0.5, 0.0]

model = 'm50n512'
simba = 's50j7k'

mmin = 8.7
mmax = 12.5

# one subplot per redshift

#fig, axs = plt.subplots(3, 2)
#axes = axs.reshape(-1)

fig = plt.figure(figsize=(15, 17))

for i in range(6):
	
	with h5py.File(data_dir+model+'_fh_qr_'+mufasa_snaps[i]+'_data.h5', 'r') as mufasa_data:
                mufasa_rad_l = mufasa_data['halflight'].value
                mufasa_central = mufasa_data['central'].value
                mufasa_smass = mufasa_data['stellar_mass'].value
                mufasa_sfr = mufasa_data['sfr'].value


	with h5py.File(data_dir+model+'_'+simba+'_'+simba_snaps[i]+'_data.h5', 'r') as simba_data:
		simba_rad_l = simba_data['halflight'].value
		simba_central = simba_data['central'].value
		simba_smass = simba_data['stellar_mass'].value
		simba_sfr = simba_data['sfr'].value
	
	simba_ssfr = 1.e9*simba_sfr/simba_smass
	simba_ssfr = np.log10(simba_ssfr+10**(-2.5+0.3*simba_z[i]))
	ssfrlim = min(simba_ssfr)+0.2

	# plot simba galaxy points
	massbin,cvecbin,ebinlo,ebinhi = pm.runningmedian(np.log10(simba_smass[simba_central]),simba_ssfr[simba_central])
        simba_x = np.log10(simba_smass[simba_rad_l>0])  
        simba_y = np.log10(simba_rad_l[simba_rad_l>0])
        simba_c = simba_ssfr - np.interp(np.log10(simba_smass),massbin,cvecbin)

	cmap = plt.get_cmap('jet_r')
	pixsize = 3*(simba_x/min(simba_x)+1)
	
	ax = fig.add_subplot(3, 2, i+1)

	im = ax.scatter(simba_x, simba_y, c=simba_c, s=pixsize, lw=0, cmap=cmap, label='Simba')
	cbar = fig.colorbar(im,ax=ax, label=r'$\Delta\log$sSFR')
	cbar.ax.tick_params(labelsize=10)
	cbar.set_clim(-1.6, )
	#cbar.ax.set_label(r'$\Delta\log$sSFR', fontsize=8) 
	
	# plot simba red and blue/ all galaxies
	if simba_z[i] <= 1.5:
        	bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(np.log10(simba_smass[simba_ssfr>ssfrlim]),np.log10(simba_rad_l[simba_ssfr>ssfrlim]))
                ax.plot(bin_cent,ymean,'--',lw=4,color='c', label='Star forming')
                
		bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(np.log10(simba_smass[simba_ssfr<ssfrlim]),np.log10(simba_rad_l[simba_ssfr<ssfrlim]))
                ax.plot(bin_cent,ymean,'--',lw=4,color='m', label='Passive')
	else:
		bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(simba_x,simba_y)
		ax.plot(bin_cent,ymean,'--',lw=4,color='b', label='Simba')
	
	# plot mufasa median data		
	mufasa_x = np.log10(mufasa_smass[mufasa_rad_l>0])
	mufasa_y = np.log10(mufasa_rad_l[mufasa_rad_l>0])
	bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(mufasa_x, mufasa_y)
	ax.plot(bin_cent,ymean,'--',lw=2,color='g', label='Mufasa')
	
	# plot redshift observational data
	plot_data(ax, simba_z[i])

	ax.annotate('z=%g'%(np.round(simba_z[i],1)), xy=(0.1, 0.9), xycoords='axes fraction',size=16,bbox=dict(boxstyle="round", fc="w"))
	ax.minorticks_on()
	ax.set_xlim(mmin,mmax)
	ax.set_ylim(-0.3,1.8-0.3*simba_z[i])
	ax.set_xlabel(r'$\log\ M_{*}$',fontsize=16)
	ax.set_ylabel(r'$\log\ R_{half,*}$' ,fontsize=16)
	ax.legend(loc='lower right', fontsize=8)
	
plt.savefig(plot_dir+'halflight_'+model+'.pdf', bbox_inches='tight', format='pdf')
plt.clf()

