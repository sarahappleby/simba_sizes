import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys
from plotting_methods import *

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 14})

softening = 0.5  * 2.8# kpc/h

data_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/data/'
plot_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/plots/'

bin_labels = ['10.0-10.5', '10.5-11.0', '>11.0']
plot_labels = [r'$10.0 < \textrm{log} (M_* / M_{\odot}) < 10.5$', r'$10.5 < \textrm{log} (M_* / M_{\odot}) < 11.0$', r'$ \textrm{log} (M_* / M_{\odot}) > 11.0$']

cmap = cm.get_cmap('plasma')
colors_b18 = [cmap(0.65), cmap(0.5), cmap(0.3)]
colors_b18 = [cmap(0.8), cmap(0.6), cmap(0.3)]

cmap = cm.get_cmap('viridis')
colors = [cmap(0.6), cmap(0.35), cmap(0.15)]
colors = [cmap(0.8), cmap(0.6), cmap(0.3)]

mass_b18=[10.0,  10.5, 11., 11.5]

fig, ax = plt.subplots(1, 2, figsize=(15, 6))

# for the star forming galaxies:
plot_belfiore(ax[0], 'sf', colors_b18, mass_b18=mass_b18, label=False)

for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+'sf_rot_ssfr_data.h5') as f:
        tukey = f['all_tukey_'+b].value
        err = f['all_err_'+b].value
    if i == 0:
        n = len(np.array(tukey))
        dr = 0.2
        factor = dr*n
        bins = np.arange(0., factor, dr)
        rplot = bins+(dr*0.5)

    ax[0].plot(rplot, tukey, color=colors[i], linestyle='-', label=plot_labels[i])
    if i < 2:
        alpha = 0.25
    elif i == 2:
        alpha = 0.15
    ax[0].fill_between(rplot, tukey - err, tukey + err, color=colors[i], alpha=alpha)
    if i == 2:
        with h5py.File(data_dir+'sf_rand_ssfr_data.h5') as f:
            tukey = f['all_tukey_'+b].value
            err = f['all_err_'+b].value
        ax[0].plot(rplot, tukey, color=colors[i], linestyle='--', label=plot_labels[i] + '; random orientation')
        ax[0].fill_between(rplot, tukey - err, tukey + err, color=colors[i], alpha=alpha)



ax[0].annotate(r'$\textbf{SF}$', xy=(0.1, 0.1), xycoords='axes fraction',size=16)
ax[0].set_xlabel(r'$R / R_{half}$', fontsize=16)
ax[0].set_ylabel(r'$\textrm{log} (\textrm{sSFR} / \textrm{yr}^{-1})$', fontsize=16)
ax[0].set_xlim(0, 2.)
ax[0].set_ylim(-13.5, -9.5)

ax[0].legend(fontsize=12, loc=4)
plot_res_limit(ax[0], softening, 'star_forming', colors)

# for the green valley galaxies:
plot_belfiore(ax[1], 'gv', colors_b18, mass_b18=mass_b18, label=True)
for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+'gv_rot_ssfr_data.h5') as f:
        tukey = f['all_tukey_'+b].value
        err = f['all_err_'+b].value
    if i == 0:
        n = len(np.array(tukey))
        dr = 0.2
        factor = dr*n
        bins = np.arange(0., factor, dr)
        rplot = bins+(dr*0.5)

    ax[1].plot(rplot, tukey, color=colors[i], linestyle='-')
    if i < 2:
        alpha = 0.25
    elif i == 2:
        alpha = 0.15
    ax[1].fill_between(rplot, tukey - err, tukey + err, color=colors[i], alpha=alpha)
    if i == 2:
        with h5py.File(data_dir+'gv_rand_ssfr_data.h5') as f:
            tukey = f['all_tukey_'+b].value
            err = f['all_err_'+b].value
        ax[1].plot(rplot, tukey, color=colors[i], linestyle='--')
        ax[1].fill_between(rplot, tukey - err, tukey + err, color=colors[i], alpha=alpha)


ax[1].annotate(r'$\textbf{GV}$', xy=(0.1, 0.9), xycoords='axes fraction',size=16)
ax[1].set_xlabel(r'$R / R_{half}$', fontsize=16)
ax[1].set_ylabel(r'$\textrm{log} (\textrm{sSFR} / \textrm{yr}^{-1})$', fontsize=16)
ax[1].set_xlim(0, 2.0)
ax[1].set_ylim(-13.5, -9.5)
ax[1].legend(fontsize=12)

plot_res_limit(ax[1], softening, 'green_valley', colors)

plt.savefig(plot_dir+'ssfr_profiles_s50j7k_belfiore.png')
#plt.clf()
