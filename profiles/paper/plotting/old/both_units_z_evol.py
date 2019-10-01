import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from plotting_methods import plot_tacchella
import sys

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 14})

cmap = cm.get_cmap('viridis')
colors = [cmap(0.8), cmap(0.6), cmap(0.3)]

cmap = cm.get_cmap('plasma')
t18_colors = [cmap(0.65), cmap(0.35)]

gals = 'all'

data_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/data/'
plot_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/plots/'
bin_labels = ['10.0-10.5', '10.5-11.0', '>11.0']
plot_labels = [r'$10.0 < \textrm{log} (M_* / M_{\odot}) < 10.5$',
                r'$10.5 < \textrm{log} (M_* / M_{\odot}) < 11.0$',
                r'$ \textrm{log} (M_* / M_{\odot}) > 11.0$']
xlabel = r'$ R / R_{half}$'

fig, ax = plt.subplots(2, 2, figsize=(12, 12))
ax = ax.flatten()


### Halflight radius units

plot_tacchella(ax[0], 'norm', t18_colors, label=True)
# plot sSFR:
for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+'snap_078_sf_rot_ssfr_data.h5') as f:
        no_gals = f[gals+'_no_gals_'+b].value
        tukey = f[gals+'_tukey_'+b].value
        small = f[gals+'_small_scale_'+b].value
        large = f[gals+'_large_scale_'+b].value
    if i == 0:
        n = len(np.array(tukey))
        dr = 0.2
        factor = dr*n
        bins = np.arange(0., factor, dr)
        rplot = bins+(dr*0.5)

    ax[0].plot(rplot, tukey, color=colors[i], linestyle='-', label=plot_labels[i]+'; '+str(int(no_gals))+' galaxies')
    if i == 1:
        ax[0].fill_between(rplot, tukey - large, tukey + large, color=colors[i], alpha=0.1)
    ax[0].fill_between(rplot, tukey - small, tukey + small, color=colors[i], alpha=0.3)

ax[0].set_ylim(-12.5, -8.0)
ax[0].legend()
ax[0].set_xlabel(xlabel)
ax[0].set_ylabel(r'$\textrm{log} (\textrm{sSFR} / \textrm{yr}^{-1})$')
ax[0].legend(fontsize=12)
ax[0].annotate('z=2', xy=(0.9, 0.9), xycoords='axes fraction',size=13,bbox=dict(boxstyle="round", fc="w"))

for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+'snap_105_sf_rot_ssfr_data.h5') as f:
        no_gals = f[gals+'_no_gals_'+b].value
        tukey = f[gals+'_tukey_'+b].value
        small = f[gals+'_small_scale_'+b].value
        large = f[gals+'_large_scale_'+b].value
    if i == 0:
        n = len(np.array(tukey))
        dr = 0.2
        factor = dr*n
        bins = np.arange(0., factor, dr)
        rplot = bins+(dr*0.5)

    ax[1].plot(rplot, tukey, color=colors[i], linestyle='-', label=str(int(no_gals))+' galaxies')
    if i == 1:
        ax[1].fill_between(rplot, tukey - large, tukey + large, color=colors[i], alpha=0.1)
    ax[1].fill_between(rplot, tukey - small, tukey + small, color=colors[i], alpha=0.3)

ax[1].set_ylim(-12.5, -8.0)
ax[1].set_xlabel(xlabel)
ax[1].set_ylabel(r'$\textrm{log} (\textrm{sSFR} / \textrm{yr}^{-1})$')
ax[1].legend(fontsize=12, loc=3)
ax[1].annotate('z=1', xy=(0.9, 0.9), xycoords='axes fraction',size=13,bbox=dict(boxstyle="round", fc="w"))

### Physical radius units

xlabel = r'$ R (\textrm{kpc})$'

plot_tacchella(ax[2], 'phys', t18_colors, label=True)
# plot sSFR:
for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+'phys_snap_078_sf_rot_ssfr_data.h5') as f:
        no_gals = f[gals+'_no_gals_'+b].value
        tukey = f[gals+'_tukey_'+b].value
        small = f[gals+'_small_scale_'+b].value
        large = f[gals+'_large_scale_'+b].value
    if i == 0:
        n = len(np.array(tukey))
        dr = 1.0
        bins = np.arange(0., n*dr, dr)
        rplot = bins+(dr*0.5)

    ax[2].plot(rplot, tukey, color=colors[i], linestyle='-')
    if i == 1:
        ax[2].fill_between(rplot, tukey - large, tukey + large, color=colors[i], alpha=0.1)
    ax[2].fill_between(rplot, tukey - small, tukey + small, color=colors[i], alpha=0.3)

ax[2].set_ylim(-12.5, -8.0)
ax[2].set_xlabel(xlabel)
ax[2].legend()
ax[2].set_ylabel(r'$\textrm{log} (\textrm{sSFR} / \textrm{yr}^{-1})$')
ax[2].annotate('z=2', xy=(0.9, 0.9), xycoords='axes fraction',size=13,bbox=dict(boxstyle="round", fc="w"))

for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+'phys_snap_105_sf_ssfr_data.h5') as f:
        no_gals = f[gals+'_no_gals_'+b].value
        tukey = f[gals+'_tukey_'+b].value
        small = f[gals+'_small_scale_'+b].value
        large = f[gals+'_large_scale_'+b].value
    if i == 0:
        n = len(np.array(tukey))
        dr = 1.0
        bins = np.arange(0., n*dr, dr)
        rplot = bins+(dr*0.5)

    ax[3].plot(rplot, tukey, color=colors[i], linestyle='-')
    if i == 1:
        ax[3].fill_between(rplot, tukey - large, tukey + large, color=colors[i], alpha=0.1)
    ax[3].fill_between(rplot, tukey - small, tukey + small, color=colors[i], alpha=0.3)

ax[3].set_ylim(-12.5, -8.0)
ax[3].set_xlabel(xlabel)
ax[3].set_ylabel(r'$\textrm{log} (\textrm{sSFR} / \textrm{yr}^{-1})$')
ax[3].annotate('z=1', xy=(0.9, 0.9), xycoords='axes fraction',size=13,bbox=dict(boxstyle="round", fc="w"))



plt.savefig(plot_dir+'both_redshift_ssfr.png')
plt.show()
plt.clf()

