import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt 
import h5py

def tukey_biweight(x, c=5.0, epsilon=1.e-11):
    median = np.nanpercentile(x, '50', axis=0)
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

def plot_belfiore(ax, sample, colors, mass_b18=[9.0,  9.5, 10.0,  10.5, 11., 11.5]):
    if sample == 'sf':
        aa = ascii.read('b18_ssfr_gradients_main_sequence.dat')
    if sample == 'gv':
        aa = ascii.read('b18_ssfr_gradients_GV.dat')
    x=aa['R']
    nmass = len(mass_b18)
    for nn in range(0, nmass-1):
        y=np.array(aa['{:.2f}'.format(mass_b18[nn])+'-'+'{:.2f}'.format(mass_b18[nn+1])+' sSFR'])
        y1=np.array(aa['{:.2f}'.format(mass_b18[nn])+'-'+'{:.2f}'.format(mass_b18[nn+1])+' sSFR']+\
            aa['{:.2f}'.format(mass_b18[nn])+'-'+'{:.2f}'.format(mass_b18[nn+1])+' error sSFR'])
        y2=np.array(aa['{:.2f}'.format(mass_b18[nn])+'-'+'{:.2f}'.format(mass_b18[nn+1])+' sSFR']-\
            aa['{:.2f}'.format(mass_b18[nn])+'-'+'{:.2f}'.format(mass_b18[nn+1])+' error sSFR'])

        #ax.plot(x,y , color=colors[nn], lw=2., label ='{:.1f}'.format(mass_b18[nn])+'-'+'{:.1f}'.format(mass_b18[nn+1]))
        ax.plot(x,y , color=colors[nn], lw=1.5)
        ax.fill_between(x, y1, y2, where=y1 >= y2, facecolor=colors[nn], alpha=0.2, edgecolor='none')

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

def plot_all(data_dirs, filename, xlabel, ylabel, xlim, ylim=None, savefile='plot.png'):

    bin_centrals = [10.0, 10.5, 11.0]
    bin_sats = [10.0, 10.5]
    labels_centrals, basic_centrals = get_labels(bin_centrals)
    labels_sats, basic_sats = get_labels(bin_sats)
    colors = ['b', 'm', 'r']
    mass_b18=[10.0,  10.5, 11.0, 11.5]

    fig, ax = plt.subplots(2, 2, figsize=(15, 15))
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
        
    if 'ssfr' in filename:
        plot_belfiore(ax[0][0], 'sf', colors, mass_b18=mass_b18)
        plot_belfiore(ax[0][1], 'gv', colors, mass_b18=mass_b18)

    plt.savefig(savefile)
    plt.clf()

