import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from plotting_methods import plot_belfiore
import sys

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 14})

cmap = cm.get_cmap('plasma')
colors_b18 = [cmap(0.85), cmap(0.65), cmap(0.5), cmap(0.3)]
colors_b18 = [cmap(0.8), cmap(0.6), cmap(0.3)]

cmap = cm.get_cmap('viridis')
colors = [cmap(0.85), cmap(0.6), cmap(0.35), cmap(0.15)]
colors = [cmap(0.8), cmap(0.6), cmap(0.3)]

#colors = ['g', 'c', 'b', 'm']
selection = 'gv'
gals = 'cen'
mass_b18 = [10.0,  10.5, 11., 11.5]

data_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/data/'
plot_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/plots/'
bin_labels = ['10.0-10.5', '10.5-11.0', '>11.0']
plot_labels = [r'$10.0 < \textrm{log} (M_* / M_{\odot}) < 10.5$', 
                r'$10.5 < \textrm{log} (M_* / M_{\odot}) < 11.0$', 
                r'$ \textrm{log} (M_* / M_{\odot}) > 11.0$']
xlabel = r'$ R / R_{half}$'

fig, ax = plt.subplots(2, 2, figsize=(13, 13))
ax = ax.flatten()

# plot SFR:
for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+selection+'_sfr_data.h5') as f:
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

    ax[0].plot(rplot, tukey, color=colors[i], linestyle='--', label=plot_labels[i]+'; '+str(int(no_gals))+' galaxies')
    if i == 1:
        ax[0].fill_between(rplot, tukey - large, tukey + large, color=colors[i], alpha=0.1)
    ax[0].fill_between(rplot, tukey - small, tukey + small, color=colors[i], alpha=0.3)

ax[0].set_ylim(-5, -1)
if gals == 'all':
    ax[0].set_xlim(0, 3.5)
else:
    ax[0].set_xlim(0., 2)

ax[0].set_xlabel(xlabel)
ax[0].set_ylabel(r'$ \textrm{log} (\Sigma_{\textrm{SFR}} / M_{\odot}\textrm{yr}^{-1} \textrm{kpc}^{-2})$')
if not gals == 'all':
    ax[0].legend(fontsize=12, loc=1)


plot_belfiore(ax[1], selection, colors_b18, mass_b18=mass_b18, label=True)
# plot sSFR:
for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+selection+'_ssfr_data.h5') as f:
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

    ax[1].plot(rplot, tukey, color=colors[i], linestyle='--')
    if i == 1:
        ax[1].fill_between(rplot, tukey - large, tukey + large, color=colors[i], alpha=0.1)
    ax[1].fill_between(rplot, tukey - small, tukey + small, color=colors[i], alpha=0.3)

if gals == 'all':
    ax[1].set_xlim(0, 3.5)
else:
    ax[1].set_xlim(0., 2)
ax[1].set_ylim(-12.5, -9.5)
ax[1].set_xlabel(xlabel)
ax[1].set_ylabel(r'$\textrm{log} (\textrm{sSFR} / \textrm{yr}^{-1})$')
if not gals == 'all':
    ax[1].legend(fontsize=12, loc=1)
else:
    ax[1].legend(fontsize=12, loc=3)

# plot HI:
for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+selection+'_h1_data.h5') as f:
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

    ax[2].plot(rplot, tukey, color=colors[i], linestyle='--', label=plot_labels[i]+'; '+str(int(no_gals))+' galaxies')
    if i == 1:
        ax[2].fill_between(rplot, tukey - large, tukey + large, color=colors[i], alpha=0.1)
    ax[2].fill_between(rplot, tukey - small, tukey + small, color=colors[i], alpha=0.3)

if not gals == 'all':
    ax[2].set_ylim(3.5, )
ax[2].set_xlim(0, 5)
ax[2].set_ylim(4.5,8.)
ax[2].set_xlabel(xlabel)
ax[2].set_ylabel(r'$ \textrm{log} (\Sigma_{HI} / M_{\odot}\textrm{kpc}^{-2})$')
if gals == 'all':
    ax[2].legend(fontsize=12, loc=3)



# plot fmol:
for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+selection+'_fmol_data.h5') as f:
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

    ax[3].plot(rplot, tukey, color=colors[i], linestyle='--')
    if i == 1:
        ax[3].fill_between(rplot, tukey - large, tukey + large, color=colors[i], alpha=0.1)
    ax[3].fill_between(rplot, tukey - small, tukey + small, color=colors[i], alpha=0.3)

ax[3].set_ylim(-2., 1.5)
ax[3].set_xlim(0, 5)
ax[3].set_xlabel(xlabel)
ax[3].set_ylabel(r'$ \textrm{log} (f_{mol})$')

plt.savefig(plot_dir+selection+'_all_radial.png')
plt.show()
plt.clf()
