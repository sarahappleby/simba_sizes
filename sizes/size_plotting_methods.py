import matplotlib.pyplot as plt
import numpy as np

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
                                ax.plot(ms_data,logRe_sf,':',color='mediumblue',lw=6,label='SDSS-SF', markersize=3)
                                alpha = -0.02
                                beta = 0.65
                                gamma = 1.58
                                M0 = 1.11e10/0.68**2
                                logRe_q = np.log10(gamma * (10**ms_data/M0)**alpha * (1+10**ms_data/M0)**(beta-alpha))
                                ax.plot(ms_data,logRe_q,':',color='r',lw=6,label='SDSS-Q', markersize=3)
                                # Kelvin+11 GAMA
                                #logRe_spiral = 0.4+0.3*(ms_data-9)
                                #plt.plot(ms_data,logRe_spiral,'-',color='k',lw=4,label='GAMA (Baldry+12)')
                if z >= 0.5 and z < 1.5: # vdWel+14 CANDELS
                                #average of z = 0.75, 1.25, table 2
                                logms = [9.25,9.75,10.25,10.75,11.25]
                                logRe = [0.4, 0.52, 0.61, 0.71, 0.86]
                                eRelo = [0.25, 0.25, 0.25, 0.22, 0.17]
                                eRehi = [0.23, 0.21, 0.09, 0.16, 0.18]
                                ax.errorbar(np.array(logms)-0.05,logRe,lw=3,yerr=[eRelo,eRehi],fmt='o',color='mediumblue', ecolor='mediumblue',label='CANDELS LTG', markersize=6)
                                logms = [9.25,9.75,10.25,10.75,11.25]
                                logRe = [0.23, 0.195, 0.165, 0.375, 0.695]
                                eRelo = [0.25, 0.33, 0.25, 0.21, 0.18]
                                eRehi = [0.2, 0.24, 0.23, 0.22, 0.20]
                                ax.errorbar(np.array(logms)+0.05,logRe,lw=3,yerr=[eRelo,eRehi],fmt='o',color='r', ecolor='r',label='CANDELS ETG', markersize=6)
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
                                ax.errorbar(np.array(logms)-0.05,logRe,lw=3,yerr=[eRelo,eRehi],fmt='o',color='mediumblue', ecolor='mediumblue',label='CANDELS LTG', markersize=6)
                                logms = [9.75, 10.25, 10.75, 11.25]
                                logRe = [0.22, -0.01, 0.14, 0.41]
                                eRelo = [0.24, 0.31, 0.26, 0.19]
                                eRehi = [0.26, 0.36, 0.38, 0.24]
                                ax.errorbar(np.array(logms)+0.05,logRe,lw=3,yerr=[eRelo,eRehi],fmt='o',color='r',ecolor='r',label='CANDELS ETG', markersize=6)
                                """ ltg + etg, table 5
                                logms = [9.75,10.25,10.75,11.25]
                                logRe = [0.39,0.44,0.47,0.62]
                                eRelo = [0.27,0.32,0.39,0.30]
                                eRehi = [0.23,0.21,0.24,0.21]
                                ax.plot(logms,logRe,'o',color='k',label='CANDELS (van der Wel+14)', markersize=6)
                                ax.errorbar(logms,logRe,lw=3,yerr=[eRelo,eRehi],color='k')
                                """

