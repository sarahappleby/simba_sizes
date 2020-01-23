import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import h5py

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
        new_cmap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
                                                                                                                cmap(np.linspace(minval, maxval, n)))
        return new_cmap

def plot_sdss_sf(ax):
    ms_data = np.linspace(9,12.0,20)
    alpha = 0.23
    beta = 0.41
    gamma = 10.72
    M0 = 4.9e11/0.68**2
    logRe_sf = np.log10(gamma * (10**ms_data/M0)**alpha * (1+10**ms_data/M0)**(beta-alpha))
    ax.plot(ms_data,logRe_sf,'-',color='dodgerblue',lw=3.5,label='SDSS-SF', markersize=3)


def plot_data(ax, z):
                if z < 0.5: # Zhang+19 SDSS
                                ms_data = np.linspace(9,12.0,20)
                                alpha = 0.23
                                beta = 0.41
                                gamma = 10.72
                                M0 = 4.9e11/0.68**2
                                logRe_sf = np.log10(gamma * (10**ms_data/M0)**alpha * (1+10**ms_data/M0)**(beta-alpha))
                                ax.plot(ms_data,logRe_sf,'-',color='dodgerblue',lw=4,label='SDSS-SF', markersize=3)
                                alpha = -0.02
                                beta = 0.65
                                gamma = 1.58
                                M0 = 1.11e10/0.68**2
                                logRe_q = np.log10(gamma * (10**ms_data/M0)**alpha * (1+10**ms_data/M0)**(beta-alpha))
                                ax.plot(ms_data,logRe_q,'-',color='r',lw=3.5,label='SDSS-Q', markersize=3)
                                # Kelvin+11 GAMA
                                #logRe_spiral = 0.4+0.3*(ms_data-9)
                                #plt.plot(ms_data,logRe_spiral,'-',color='k',lw=4,label='GAMA (Baldry+12)')
                if z >= 0.5 and z < 1.5: # vdWel+14 CANDELS
                                #average of z = 0.75, 1.25, table 2
                                logms = [9.25,9.75,10.25,10.75,11.25]
                                logRe = [0.4, 0.52, 0.61, 0.71, 0.86]
                                eRelo = [0.25, 0.25, 0.25, 0.22, 0.17]
                                eRehi = [0.23, 0.21, 0.09, 0.16, 0.18]
                                ax.errorbar(np.array(logms)-0.05,logRe,lw=3,yerr=[eRelo,eRehi],fmt='s',color='dodgerblue', ecolor='dodgerblue',label='CANDELS LTG', markersize=6)
                                logms = [9.25,9.75,10.25,10.75,11.25]
                                logRe = [0.23, 0.195, 0.165, 0.375, 0.695]
                                eRelo = [0.25, 0.33, 0.25, 0.21, 0.18]
                                eRehi = [0.2, 0.24, 0.23, 0.22, 0.20]
                                ax.errorbar(np.array(logms)+0.05,logRe,lw=3,yerr=[eRelo,eRehi],fmt='s',color='r', ecolor='r',label='CANDELS ETG', markersize=6)
                if z >= 1.5 and z < 2.5:
                                # Allen+17 CANDELS+ZFOURGE, table 1
                                ms_data = np.linspace(9,11.5,5)
                                re_data = 0.2*(ms_data-10)+0.4
                                ax.plot(ms_data,re_data,'-',color='k',lw=2.5,label='CANDELS+ZFOURGE')
                                # vdWel+14 CANDELS, table 2
                                # average of z = 1.75, 2.25, table 2
                                logms = [9.75, 10.25, 10.75, 11.25]
                                logRe = [0.39, 0.48, 0.57, 0.67]
                                eRelo = [0.26, 0.26, 0.27, 0.21]
                                eRehi = [0.22, 0.2, 0.18, 0.19]
                                ax.errorbar(np.array(logms)-0.05,logRe,lw=3,yerr=[eRelo,eRehi],fmt='s',color='dodgerblue', ecolor='dodgerblue',label='CANDELS LTG', markersize=6)
                                logms = [9.75, 10.25, 10.75, 11.25]
                                logRe = [0.22, -0.01, 0.14, 0.41]
                                eRelo = [0.24, 0.31, 0.26, 0.19]
                                eRehi = [0.26, 0.36, 0.38, 0.24]
                                ax.errorbar(np.array(logms)+0.05,logRe,lw=3,yerr=[eRelo,eRehi],fmt='s',color='r',ecolor='r',label='CANDELS ETG', markersize=6)
                                """ ltg + etg, table 5
                                logms = [9.75,10.25,10.75,11.25]
                                logRe = [0.39,0.44,0.47,0.62]
                                eRelo = [0.27,0.32,0.39,0.30]
                                eRehi = [0.23,0.21,0.24,0.21]
                                ax.plot(logms,logRe,'o',color='k',label='CANDELS (van der Wel+14)', markersize=6)
                                ax.errorbar(logms,logRe,lw=3,yerr=[eRelo,eRehi],color='k')
                                """

def plot_high_res(ax, snap, mtype):

    high_res_file = '/home/sapple/simba_sizes/sizes/data/m25n512_pyloser_medians.h5'

    with h5py.File(high_res_file, 'r') as hrf:

        # first do star forming galaxies
        bin_cent = hrf['snap_'+str(snap)][mtype+'_sf_bin_cent']
        median = hrf['snap_'+str(snap)][mtype+'_sf_median']
        ysighi = hrf['snap_'+str(snap)][mtype+'_sf_ysighi']
        ysiglo = hrf['snap_'+str(snap)][mtype+'_sf_ysiglo']

        (_, caps, _) = ax.errorbar(bin_cent, median, yerr=[ysiglo, ysighi], fmt='o', ls='-', lw=2.5, color='blue', 
                capsize=7, label='Star forming (high res)')

        for cap in caps:
            cap.set_markeredgewidth(1)

        # next quenched galaxies

        bin_cent = hrf['snap_'+str(snap)][mtype+'_q_bin_cent']
        median = hrf['snap_'+str(snap)][mtype+'_q_median']
        ysighi = hrf['snap_'+str(snap)][mtype+'_q_ysighi']
        ysiglo = hrf['snap_'+str(snap)][mtype+'_q_ysiglo']

        (_, caps, _) = ax.errorbar(bin_cent, median, yerr=[ysiglo, ysighi], fmt='o', ls='-', lw=2.5, color='m', 
                capsize=7, label='Passive (high res)')

        for cap in caps:
            cap.set_markeredgewidth(1)

def do_legend_order(ax, z):
    if z == 0.:
        labels_use = ['SDSS-Q','SDSS-SF','Star forming','Passive','Star forming (high res)','Passive (high res)']
    elif z== 1.:
        labels_use = ['CANDELS LTG','CANDELS ETG','Star forming','Passive', 'Star forming (high res)','Passive (high res)']
    elif z==2.:
        labels_use = ['CANDELS+ZFOURGE','CANDELS LTG','CANDELS ETG','Star forming','Passive','Star forming (high res)','Passive (high res)']

    handles, labels = ax.get_legend_handles_labels()
    index = [labels.index(l) for l in labels_use]
    handles = list(np.array(handles)[index])
    
    ax.legend(handles, labels_use, fontsize=11)
