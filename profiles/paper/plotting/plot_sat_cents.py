import h5py
import numpy as np
import matplotlib.pyplot as plt
from plotting_methods import plot_belfiore
import sys

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 12})

colors = ['g', 'c', 'b', 'm', 'r']
selection = 'sf'
colors_b18 = ['lightgreen', 'c', 'deepskyblue', 'magenta']
mass_b18 = [9.5, 10.0,  10.5, 11., 11.5]

data_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/data/'
plot_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/plots/'
bin_labels = ['9.5-10.0', '10.0-10.5', '10.5-11.0', '>11.0']
plot_labels = [r'$9.5 < \textrm{log} (M_* / M_{\odot}) < 10.0$',
                r'$10.0 < \textrm{log} (M_* / M_{\odot}) < 10.5$',
                r'$10.5 < \textrm{log} (M_* / M_{\odot}) < 11.0$',
                r'$ \textrm{log} (M_* / M_{\odot}) > 11.0$']
xlabel = r'$ R / R_{half}$'

fig, ax = plt.subplots(1, 2, figsize=(15, 10))
ax = ax.flatten()

# plot SFR:
for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+'sf_ssfr_data.h5') as f:
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

    ax[0].plot(rplot, tukey, color=colors[i], marker='.', markersize=4, linestyle=':', label=plot_labels[i]+'; '+str(int(no_gals))+' galaxies')
    if i == 1:
        ax[0].fill_between(rplot, tukey - large, tukey + large, color=colors[i], alpha=0.1)
    ax[0].fill_between(rplot, tukey - small, tukey + small, color=colors[i], alpha=0.3)

    with h5py.File(data_dir+'sf_ssfr_data.h5') as f:
        no_gals = f['sat_no_gals_'+b].value
        tukey = f['sat_tukey_'+b].value
        small = f['sat_small_scale_'+b].value
        large = f['sat_large_scale_'+b].value
    
    #for satellites, join last two mass bins


    if i == 3:


    ax[0].plot(rplot, tukey, color=colors[i], marker='.', markersize=4, linestyle='--', label=plot_labels[i]+'; '+str(int(no_gals))+' galaxies')
    ax[0].fill_between(rplot, tukey - small, tukey + small, color=colors[i], alpha=0.3)
ax[0].set_ylim(-12.5, -9.0)
ax[0].set_xlim(0, 3.5)
ax[0].set_xlabel(xlabel)
ax[0].set_ylabel(r'$\textrm{log} (\textrm{sSFR} / \textrm{yr}^{-1})$')
ax[0].legend()

for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+'sf_h1_data.h5') as f:
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

    ax[1].plot(rplot, tukey, color=colors[i], marker='.', markersize=4, linestyle=':', label=plot_labels[i]+'; '+str(int(no_gals))+' galaxies')
    if i == 1:
        ax[1].fill_between(rplot, tukey - large, tukey + large, color=colors[i], alpha=0.1)
    ax[1].fill_between(rplot, tukey - small, tukey + small, color=colors[i], alpha=0.3)

    with h5py.File(data_dir+'sf_h1_data.h5') as f:
        no_gals = f['sat_no_gals_'+b].value
        tukey = f['sat_tukey_'+b].value
        small = f['sat_small_scale_'+b].value
        large = f['sat_large_scale_'+b].value
    if i == 0:
        n = len(np.array(tukey))
        dr = 0.2
        factor = dr*n
        bins = np.arange(0., factor, dr)
        rplot = bins+(dr*0.5)

    ax[1].plot(rplot, tukey, color=colors[i], marker='.', markersize=4, linestyle='--', label=plot_labels[i]+'; '+str(int(no_gals))+' galaxies')
    ax[1].fill_between(rplot, tukey - small, tukey + small, color=colors[i], alpha=0.3)

ax[1].set_xlim(0, 5.)
ax[1].set_xlabel(xlabel)
ax[1].set_ylabel(r'$\textrm{log} (\Sigma_{HI} / M_{\odot}\textrm{kpc}^{-2})$')
ax[1].legend()

plt.show()

