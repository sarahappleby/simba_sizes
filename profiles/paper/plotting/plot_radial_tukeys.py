import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from plotting_methods import plot_belfiore
import sys

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 12})

cmap = cm.get_cmap('magma')
colors_b18 = [cmap(0.9), cmap(0.65), cmap(0.5), cmap(0.3)]

cmap = cm.get_cmap('viridis')
colors = [cmap(0.85), cmap(0.6), cmap(0.35), cmap(0.15)]

#colors = ['g', 'c', 'b', 'm']
selection = 'sf'
mass_b18 = [9.5, 10.0,  10.5, 11., 11.5]

data_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/data/'
plot_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/plots/'
bin_labels = ['9.5-10.0', '10.0-10.5', '10.5-11.0', '>11.0']
plot_labels = [r'$9.5 < \textrm{log} (M_* / M_{\odot}) < 10.0$',
                r'$10.0 < \textrm{log} (M_* / M_{\odot}) < 10.5$', 
                r'$10.5 < \textrm{log} (M_* / M_{\odot}) < 11.0$', 
                r'$ \textrm{log} (M_* / M_{\odot}) > 11.0$']
xlabel = r'$ R / R_{half}$'

fig, ax = plt.subplots(2, 2, figsize=(15, 15))
ax = ax.flatten()

# plot SFR:
for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+'sf_sfr_data.h5') as f:
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

    ax[0].plot(rplot, tukey, color=colors[i], marker='.', markersize=5, linestyle='-', label=plot_labels[i]+'; '+str(int(no_gals))+' galaxies')
    if i == 1:
        ax[0].fill_between(rplot, tukey - large, tukey + large, color=colors[i], alpha=0.1)
    ax[0].fill_between(rplot, tukey - small, tukey + small, color=colors[i], alpha=0.3)

ax[0].set_ylim(-5, )
ax[0].set_xlim(0, 3.5)
ax[0].set_xlabel(xlabel)
ax[0].set_ylabel(r'$ \textrm{log} (\Sigma_{\textrm{SFR}} / M_{\odot}\textrm{yr}^{-1} \textrm{kpc}^{-2})$')
ax[0].legend()

plot_belfiore(ax[1], 'sf', colors_b18, mass_b18=mass_b18, label=True)
# plot sSFR:
for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+'sf_ssfr_data.h5') as f:
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

    ax[1].plot(rplot, tukey, color=colors[i], marker='.', markersize=5, linestyle='-')
    if i == 1:
        ax[1].fill_between(rplot, tukey - large, tukey + large, color=colors[i], alpha=0.1)
    ax[1].fill_between(rplot, tukey - small, tukey + small, color=colors[i], alpha=0.3)

ax[1].set_xlim(0, 3.5)
ax[1].set_ylim(-12.5, -9.0)
ax[1].set_xlabel(xlabel)
ax[1].set_ylabel(r'$\textrm{log} (\textrm{sSFR} / \textrm{yr}^{-1})$')
ax[1].legend()

# plot HI:
for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+'sf_h1_data.h5') as f:
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

    ax[2].plot(rplot, tukey, color=colors[i], marker='.', markersize=5, linestyle='-')
    if i == 1:
        ax[2].fill_between(rplot, tukey - large, tukey + large, color=colors[i], alpha=0.1)
    ax[2].fill_between(rplot, tukey - small, tukey + small, color=colors[i], alpha=0.3)

ax[2].set_xlim(0, 5)
ax[2].set_xlabel(xlabel)
ax[2].set_ylabel(r'$ \textrm{log} (\Sigma_{HI} / M_{\odot}\textrm{kpc}^{-2})$')

# plot fmol:
for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+'sf_fmol_data.h5') as f:
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

    ax[3].plot(rplot, tukey, color=colors[i], marker='.', markersize=5, linestyle='-')
    if i == 1:
        ax[3].fill_between(rplot, tukey - large, tukey + large, color=colors[i], alpha=0.1)
    ax[3].fill_between(rplot, tukey - small, tukey + small, color=colors[i], alpha=0.3)

ax[3].set_ylim(-1.5, )
ax[3].set_xlim(0, 5)
ax[3].set_xlabel(xlabel)
ax[3].set_ylabel(r'$ \textrm{log} (f_{mol})$')

plt.savefig(plot_dir+selection+'_all_radial.png')
plt.show()
plt.clf()
