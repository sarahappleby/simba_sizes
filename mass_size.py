import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

import os

import caesar
import numpy as np

def plot_data(z):
        if z < 0.5: # Zhang+17 SDSS
                ms_data = np.linspace(9,12.0,20)
                alpha = 0.24
                beta = 1.33
                gamma = 10.17
                M0 = 6.49e11/0.68**2
                logRe_spiral = np.log10(gamma * (10**ms_data/M0)**alpha * (1+10**ms_data/M0)**(beta-alpha))
                plt.plot(ms_data,logRe_spiral,':',color='b',lw=4,label='SDSS-LTG (Zhang+17)')
                alpha = 0.17
                beta = 0.58
                gamma = 2.24
                M0 = 2.11e10/0.68**2
                logRe_etg = np.log10(gamma * (10**ms_data/M0)**alpha * (1+10**ms_data/M0)**(beta-alpha))
                plt.plot(ms_data,logRe_etg,':',color='r',lw=4,label='SDSS-ETG (Zhang+17)')
                # Kelvin+11 GAMA
                #logRe_spiral = 0.4+0.3*(ms_data-9)
                #plt.plot(ms_data,logRe_spiral,'-',color='k',lw=4,label='GAMA (Baldry+12)')
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
                plt.plot(ms_data,re_data,'-',color='k',lw=4,label='CANDELS+ZFOURGE (Alcorn+16)')
                # vdWel+14 CANDELS
                logms = [9.75,10.25,10.75,11.25]
                logRe = [0.39,0.44,0.47,0.62]
                eRelo = [0.27,0.32,0.39,0.30]
                eRehi = [0.23,0.21,0.24,0.21]
                plt.plot(logms,logRe,'o',color='k',ms=6,label='CANDELS (van der Wel+14)')
                plt.errorbar(logms,logRe,lw=2,yerr=[eRelo,eRehi],color='k')


data_dir = '/home/sapple/simba_sizes/data/'
plot_dir = '/home/sapple/simba_sizes/plots/'

mufasa_snaps = ['070', '085', '095', '105', '125', '126', '135']
simba_snaps = ['062', '078', '090', '105', '125', '145', '151']

mufasa_z = [3.0, 2.0, 1.5, 1.0, 0.5, 0.2, 0.0]
simba_z = [3.0, 2.0, 1.5, 1.0, 0.5, 0.1, 0.0]



# one subplot per redshift

fig, axs = plt.subplots(3, 2)
axes = axs.reshape(-1)

for i in range(len(axes)):
	simba_rad_l = np.load(data_dir+'halflight_m50n512_s50j7k_'+simba_snaps[i]+'.npy')
	simba_smass = np.load(data_dir+'s_mass_m50n512_s50j7k_'+simba_snaps[i]+'.npy')
	if not os.path.isfile(data_dir+'ssfr')

	mufasa_rad_l = np.load(data_dir+'halflight_m50n512_fh_qr_'+mufasa_snaps[i]+'.npy')
        mufasa_smass = np.load(data_dir+'s_mass_m50n512_fh_qr_'+mufasa_snaps[i]+'.npy')

	simba_x = np.log10(simba_smass[simba_rad_l>0])	
	simba_y = np.log10(simba_rad_l[simba_rad_l>0])
	
	
