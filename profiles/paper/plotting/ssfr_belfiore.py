import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys
from plotting_methods import *

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 14})

data_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/data/'
plot_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/plots/'

masks = [2, 3, 4]
bin_labels = ['10.0-10.5', '10.5-11.0', '>11.0']
plot_labels = [r'$10.0 < \textrm{log} (M_* / M_{\odot}) < 10.5$', r'$10.5 < \textrm{log} (M_* / M_{\odot}) < 11.0$', r'$ \textrm{log} (M_* / M_{\odot}) > 11.0$']

cmap = cm.get_cmap('plasma')
colors_b18 = [cmap(0.65), cmap(0.5), cmap(0.3)]
colors_b18 = [cmap(0.8), cmap(0.6), cmap(0.3)]

cmap = cm.get_cmap('viridis')
colors = [cmap(0.6), cmap(0.35), cmap(0.15)]
colors = [cmap(0.8), cmap(0.6), cmap(0.3)]

angle = 'rot'
mass_b18=[10.0,  10.5, 11., 11.5]

fig, ax = plt.subplots(1, 2, figsize=(15, 6))

# for the star forming galaxies:
plot_belfiore(ax[0], 'sf', colors_b18, mass_b18=mass_b18, label=False)

for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+'sf_'+angle+'_ssfr_data.h5') as f:
        no_gals = f['all_no_gals_'+b].value
        tukey = f['all_tukey_'+b].value
        small = f['all_small_scale_'+b].value
        large = f['all_large_scale_'+b].value
    if i == 0:
        n = len(np.array(tukey))
        dr = 0.2
        factor = dr*n
        bins = np.arange(0., factor, dr)
        rplot = bins+(dr*0.5)

    ax[0].plot(rplot, tukey, color=colors[i], linestyle='-', label=plot_labels[i] +', '+str(int(no_gals))+' galaxies')
    if i == 1:
        ax[0].fill_between(rplot, tukey - large, tukey + large, color=colors[i], alpha=0.1, edgecolor='face')
    ax[0].fill_between(rplot, tukey - small, tukey + small, color=colors[i], alpha=0.3, edgecolor='face')

ax[0].set_xlabel(r'$R / R_{half}$', fontsize=16)
ax[0].set_ylabel(r'$\textrm{log} (\textrm{sSFR} / \textrm{yr}^{-1})$', fontsize=16)
ax[0].set_xlim(0, 3.5)
ax[0].set_ylim(-12.5, -9.5)
ax[0].legend(fontsize=12)

# for the green valley galaxies:
plot_belfiore(ax[1], 'gv', colors_b18, mass_b18=mass_b18, label=True)
for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+'gv_'+angle+'_ssfr_data.h5') as f:
        no_gals = f['cen_no_gals_'+b].value
        tukey = f['cen_tukey_'+b].value
        small = f['cen_small_scale_'+b].value
        large = f['cen_large_scale_'+b].value
    if i == 0:
        n = len(np.array(tukey))
        dr = 0.2
        factor = dr*n
        bins = np.arange(0., factor, dr)
        rplot = bins+(dr*0.5)

    ax[1].plot(rplot, tukey, color=colors[i], linestyle='-')
    if i == 1:
        ax[1].fill_between(rplot, tukey - large, tukey + large, color=colors[i], alpha=0.1, edgecolor='face')
    ax[1].fill_between(rplot, tukey - small, tukey + small, color=colors[i], alpha=0.3, edgecolor='face')

ax[1].set_xlabel(r'$R / R_{half}$', fontsize=16)
ax[1].set_ylabel(r'$\textrm{log} (\textrm{sSFR} / \textrm{yr}^{-1})$', fontsize=16)
ax[1].set_xlim(0, 2.0)
ax[1].set_ylim(-12.5, -9.5)
ax[1].legend(fontsize=12)

plt.savefig(plot_dir+'ssfr_profiles_s50j7k_belfiore.png')
plt.show()
plt.clf()
