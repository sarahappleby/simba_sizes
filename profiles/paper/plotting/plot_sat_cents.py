import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import cm
from plotting_methods import plot_belfiore
import sys

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 12})

cmap = cm.get_cmap('viridis')
colors = [cmap(0.85), cmap(0.6), cmap(0.35), cmap(0.15)]

plot_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/plots/'
selection = 'sf'

data_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/data/'
plot_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/plots/'
bin_labels = ['9.5-10.0', '10.0-10.5', '>10.5']
plot_labels = [r'$9.5 < \textrm{log} (M_* / M_{\odot}) < 10.0$',
                r'$10.0 < \textrm{log} (M_* / M_{\odot}) < 10.5$',
                r'$ \textrm{log} (M_* / M_{\odot}) > 10.5$']
xlabel = r'$ R / R_{half}$'

lines = []
for c in colors:
    lines.append(Line2D([0,1],[0,1],linestyle='-', color=c))

line_cen = Line2D([0,1],[0,1],linestyle=':', color='k')
line_sat = Line2D([0,1],[0,1],linestyle='--', color='k')

fig, ax = plt.subplots(1, 2, figsize=(12, 5))
ax = ax.flatten()

leg1 = ax[0].legend(lines, plot_labels, loc=1)
leg2 = ax[0].legend([line_cen, line_sat],['centrals', 'satellites'], loc=3)
ax[0].add_artist(leg1)

leg1 = ax[1].legend(lines, plot_labels, loc=1)
leg2 = ax[1].legend([line_cen, line_sat],['centrals', 'satellites'], loc=3)
ax[1].add_artist(leg1)

# plot sSFR:
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

    ax[0].plot(rplot, tukey, color=colors[i], marker='.', markersize=4, linestyle=':')
    if i == 1:
        ax[0].fill_between(rplot, tukey - large, tukey + large, color=colors[i], alpha=0.1, edgecolor='none')
    ax[0].fill_between(rplot, tukey - small, tukey + small, color=colors[i], alpha=0.3, edgecolor='none')

    with h5py.File(data_dir+'sf_ssfr_data.h5') as f:
        no_gals = f['sat_no_gals_'+b].value
        tukey = f['sat_tukey_'+b].value
        small = f['sat_small_scale_'+b].value
        large = f['sat_large_scale_'+b].value

    ax[0].plot(rplot, tukey, color=colors[i], marker='.', markersize=4, linestyle='--')
    ax[0].fill_between(rplot, tukey - small, tukey + small, color=colors[i], alpha=0.3, edgecolor='none')
ax[0].set_ylim(-12.5, -9.0)
ax[0].set_xlim(0, 3.5)
ax[0].set_xlabel(xlabel)
ax[0].set_ylabel(r'$\textrm{log} (\textrm{sSFR} / \textrm{yr}^{-1})$')


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

    ax[1].plot(rplot, tukey, color=colors[i], marker='.', markersize=4, linestyle=':')
    if i == 1:
        ax[1].fill_between(rplot, tukey - large, tukey + large, color=colors[i], alpha=0.1, edgecolor='none')
    ax[1].fill_between(rplot, tukey - small, tukey + small, color=colors[i], alpha=0.3, edgecolor='none')

    with h5py.File(data_dir+'sf_h1_data.h5') as f:
        no_gals = f['sat_no_gals_'+b].value
        tukey = f['sat_tukey_'+b].value
        small = f['sat_small_scale_'+b].value
        large = f['sat_large_scale_'+b].value
    ax[1].plot(rplot, tukey, color=colors[i], marker='.', markersize=4, linestyle='--')
    ax[1].fill_between(rplot, tukey - small, tukey + small, color=colors[i], alpha=0.3, edgecolor='none')

ax[1].set_xlim(0, 5.)
ax[1].set_xlabel(xlabel)
ax[1].set_ylabel(r'$\textrm{log} (\Sigma_{HI} / M_{\odot}\textrm{kpc}^{-2})$')


plt.savefig(plot_dir+'sf_sats_cents.png')
plt.show()
plt.clf()
