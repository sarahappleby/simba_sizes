import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt 
import matplotlib.cm as cm
import h5py

def tukey_biweight(x, c=5.0, epsilon=1.e-11):
    median = np.nanpercentile(x, 50, axis=0)
    mad = np.nanmedian(np.abs(x - median), axis=0)
    u = (x - median)/(c*mad + epsilon)
    weights = np.zeros(u.shape)
    mask = np.abs(u) < 1.
    weights[mask] = (1. - u[mask]**2)**2
    tukey = np.nansum(x*weights, axis=0) / np.sum(weights, axis=0)
    n = x.shape[0]
    num = np.sqrt(np.nansum(((x - tukey)**2. * (1. - u**2.)**4.)*mask, axis=0))
    den = np.abs(np.nansum(((1. - u )*(1 - 5.*(u**2.)))*mask, axis=0))
    scale = np.sqrt(n)*num / den 
    t = 1.96 # 97.5% of t distribution with df = max(0.7(n-1), 1)
    return tukey, scale

def plot_kennicutt(ax, color='k', labels=False):
    def kennicutt_schmidt(sigma_gas):
        # sigma_gas in Msun/pc**2
        return (1.4 *sigma_gas) -3.6
    sigma_gas = np.arange(1, 3.2, 0.2)
    sigma_sfr = kennicutt_schmidt(sigma_gas)
    sigma_sfr = np.log10((10.**sigma_sfr) / 1.75)
    
    #low_x = np.array([0.17, 0.8]) # points taken from Kennicutt and Evans 2012, fig 12
    #low_y = np.array([-4.5, -2.46])
    #low_y = np.log10((10.**low_y) / 1.75) # correct by 1.7 to convert from Salpeter to Chabrier IMF

    if not labels:
        #ax.plot(low_x, low_y, ls='-', c=color, lw=2)
        ax.plot(sigma_gas, sigma_sfr, ls='--', c=color, lw=2)
    else:
        ax.plot(sigma_gas, sigma_sfr, ls='--', c=color, label=labels[0], lw=2)
        #ax.plot(low_x, low_y, ls='-', c=color, lw=2, label=labels[1])
    return

def plot_tacchella(ax, radial, colors, label=False):
    if radial == 'phys':
        aa = ascii.read('t18_sSFR_physR.dat')
    elif radial == 'norm':
        aa = ascii.read('t18_sSFR_normR.dat')
    x = aa['radius']
    lowM_50 = np.log10(np.array(aa['sSFR_lowM_50']) /1.e9)
    lowM_16 = np.log10(np.array(aa['sSFR_lowM_16']) /1.e9)
    lowM_84 = np.log10(np.array(aa['sSFR_lowM_84']) /1.e9)
    highM_50 = np.log10(np.array(aa['sSFR_highM_50']) /1.e9)
    highM_16 = np.log10(np.array(aa['sSFR_highM_16']) /1.e9)
    highM_84 = np.log10(np.array(aa['sSFR_highM_84']) /1.e9)
    if not label:
        ax.plot(x,lowM_50 , color=colors[0], lw=1.5)
        ax.plot(x,highM_50, color=colors[1], lw=1.5)
    elif label:
        ax.plot(x,lowM_50 , color=colors[0], lw=1.5, label=r'$\textbf{T18};\ 10.0 < \textrm{log} (M_* / M_{\odot}) < 11.0$')
        ax.plot(x,highM_50 , color=colors[1], lw=1.5, label=r'$\textbf{T18};\ \textrm{log} (M_* / M_{\odot}) > 11.0$')
    ax.fill_between(x, lowM_84, lowM_16, where=lowM_84 >= lowM_16, facecolor=colors[0], alpha=0.25, edgecolor='none')
    ax.fill_between(x, highM_84, highM_16, where=highM_84 >= highM_16, facecolor=colors[1], alpha=0.25, edgecolor='none')


def plot_belfiore(ax, sample, colors, mass_b18=[9.0,  9.5, 10.0,  10.5, 11., 11.5], label=False):
    if sample == 'sf':
        aa = ascii.read('b18_ssfr_gradients_main_sequence.dat')
    elif sample == 'gv':
        aa = ascii.read('b18_ssfr_gradients_GV.dat')
    x=aa['R']
    nmass = len(mass_b18)
    for nn in range(0, nmass-1):
        y=np.array(aa['{:.2f}'.format(mass_b18[nn])+'-'+'{:.2f}'.format(mass_b18[nn+1])+' sSFR'])
        y1=np.array(aa['{:.2f}'.format(mass_b18[nn])+'-'+'{:.2f}'.format(mass_b18[nn+1])+' sSFR']+\
            aa['{:.2f}'.format(mass_b18[nn])+'-'+'{:.2f}'.format(mass_b18[nn+1])+' error sSFR'])
        y2=np.array(aa['{:.2f}'.format(mass_b18[nn])+'-'+'{:.2f}'.format(mass_b18[nn+1])+' sSFR']-\
            aa['{:.2f}'.format(mass_b18[nn])+'-'+'{:.2f}'.format(mass_b18[nn+1])+' error sSFR'])

        if not label:
            ax.plot(x,y , color=colors[nn], lw=1.5, linestyle='--')
        elif label:
            ax.plot(x,y , color=colors[nn], lw=1.5, linestyle='--', label=r'$\textbf{B18};\ $'+ '{:.1f}'.format(mass_b18[nn])+r'$\ < \textrm{log} (M_* / M_{\odot}) <\ $'+'{:.1f}'.format(mass_b18[nn+1]))
        ax.fill_between(x, y1, y2, where=y1 >= y2, facecolor=colors[nn], alpha=0.3, edgecolor='none')
    return

def plot_res_limit(ax, softening, gal_type, colors, choose_mass_bin=None):
    median_data = {'star_forming': {'order': ['low', 'int', 'high'],
                                    'sizes': [3.4, 4.7, 5.5], 
                                    'stellar_mass': [10.2, 10.7, 11.1]},
                   'green_valley': {'order': ['low', 'int', 'high'],
                                    'sizes': [4.0, 4.8, 5.7],
                                    'stellar_mass': [10.2, 10.7, 11.1]}
                    }
    if choose_mass_bin:
        ind = [median_data['star_forming']['order'].index(choose_mass_bin)] 
    else:
        ind = [0, 1, 2]

    for i in ind:
            res = softening / median_data[gal_type]['sizes'][i]
            ax.axvline(res, c=colors[i], ls='--', lw=1.)
    return 



def get_labels(bins):
    text = '\\textrm{log} (M_* / M_{\\odot})'
    plot_labels = []
    h5py_labels = []
    for i in range(len(bins)):
        if i < len(bins) -1:
            plot_labels.append(r'$' + str(bins[i]) + ' < ' + text + ' < ' + str(bins[i+1]) + '$')
            h5py_labels.append(str(bins[i]) + '-'+str(bins[i+1]))
        else:
            plot_labels.append(r'$' + text + ' > ' + str(bins[i]) + '$')
            h5py_labels.append('>'+str(bins[i]))
    return plot_labels, h5py_labels

def plot_all(data_dirs, filename, xlabel, ylabel, xlim, ylim=None, savefile='plot.png', h1_coldens=False, h2_coldens=False):

    bin_centrals = [10.0, 10.5, 11.0]
    bin_sats = [10.0, 10.5]
    labels_centrals, basic_centrals = get_labels(bin_centrals)
    labels_sats, basic_sats = get_labels(bin_sats)
    colors = ['b', 'm', 'r']
    mass_b18=[10.0,  10.5, 11.0, 11.5]

    fig, ax = plt.subplots(2, 2, figsize=(17, 15))
    axes = ax.flatten()

    for i in range(len(axes)):
        if i <= 1:
            labels = labels_centrals
            basic = basic_centrals
        else:
            labels = labels_sats
            basic = basic_sats
        no_gals = []
        tukey = []
        large = []
        small = []
        with h5py.File(data_dirs[i]+filename, 'r') as f:
            for b in basic:
                no_gals.append(f['no_gals_'+b].value)
                tukey.append(f['tukey_'+b].value)
                large.append(f['large_scale_'+b].value)
                small.append(f['small_scale_'+b].value)

        if i == 0:
            n = np.array(tukey).shape[1]
            dr = 0.2
            factor = dr*n
            bins = np.arange(0., factor, dr)

        for m in range(len(labels)):
            axes[i].plot(bins+(dr*0.5), tukey[m], color=colors[m], marker='.', markersize=4, linestyle='--', label=labels[m] +'; '+str(int(no_gals[m]))+' galaxies')
            if m == 0:
                axes[i].fill_between(bins+(dr*0.5), tukey[m] - large[m], tukey[m] + large[m], color=colors[m], alpha=0.1)
            axes[i].fill_between(bins+(dr*0.5), tukey[m] - small[m], tukey[m] + small[m], color=colors[m], alpha=0.3)

        axes[i].set_xlim(xlim[0], xlim[1])
        if ylim:
            axes[i].set_ylim(ylim[0], ylim[1])
        axes[i].set_xlabel(xlabel, fontsize=16)
        axes[i].set_ylabel(ylabel, fontsize=16)
        axes[i].legend()

        if h1_coldens:
            lower = 3.
            higher = 9.0
            axes[i].set_ylim(lower, higher)
            axes[i].axhline(6., linestyle='--', c='k')
            ax2 = axes[i].twinx()
            convert = 1.24e14
            ax2.set_ylim(np.log10(convert*(10**lower)), np.log10(convert*(10**higher)))
            ax2.set_ylabel(r'$ \textrm{log} (N_{HI} / cm^{-2})$')
        elif h2_coldens:
            lower = 3.
            higher = 9.0
            axes[i].set_ylim(lower, higher)
            ax2 = axes[i].twinx()
            convert = 1.24e14
            ax2.set_ylim(np.log10(convert*(10**lower)), np.log10(convert*(10**higher)))
            ax2.set_ylabel(r'$ \textrm{log} (N_{H_{2}} / cm^{-2})$')
        
    if 'ssfr' in filename:
        cmap = cm.get_cmap('viridis') 
        colors_b18 = [cmap(0.3), cmap(0.6), cmap(0.9)]
        plot_belfiore(ax[0][0], 'sf', colors_b18, mass_b18=mass_b18)
        plot_belfiore(ax[0][1], 'gv', colors_b18, mass_b18=mass_b18)

    plt.savefig(savefile)
    plt.clf()

