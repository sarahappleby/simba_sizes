import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import cm
from plotting_methods import plot_belfiore
import sys

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 14})

cmap = cm.get_cmap('viridis')
colors = [cmap(0.85), cmap(0.6), cmap(0.35), cmap(0.15)]
colors = [cmap(0.6), cmap(0.3)]

plot_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/plots/'
selection = 'sf'

data_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/data/'
plot_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/plots/'
bin_labels = ['10.0-10.6', '>10.6']
plot_labels = [r'$10.0 < \textrm{log} (M_* / M_{\odot}) < 10.6$',
                r'$ \textrm{log} (M_* / M_{\odot}) > 10.6$']
xlabel = r'$ R / R_{half}$'

lines = []
for c in colors:
    lines.append(Line2D([0,1],[0,1],linestyle='-', color=c))

line_cen = Line2D([0,1],[0,1],linestyle='-', color='k')
line_sat = Line2D([0,1],[0,1],linestyle=':', color='k')

fig, ax = plt.subplots(1, 2, figsize=(12, 5))
ax = ax.flatten()

leg1 = ax[0].legend(lines, plot_labels, loc=1, fontsize=12)
leg2 = ax[0].legend([line_cen, line_sat],['centrals', 'satellites'], loc=3, fontsize=12)
ax[0].add_artist(leg1)

leg1 = ax[1].legend(lines, plot_labels, loc=1, fontsize=12)
leg2 = ax[1].legend([line_cen, line_sat],['centrals', 'satellites'], loc=3, fontsize=12)
ax[1].add_artist(leg1)

# plot sSFR:
for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+'s50j7k_star_forming_rotated_faceon_sat_cent_ratio_data.h5') as f:
        ssfr_ratio = f['ssfr_ratio_'+b].value
        ssfr_err = np.abs(f['ssfr_err_'+b].value)
    if i == 0:
        n = len(np.array(ssfr_ratio))
        dr = 0.2
        factor = dr*n
        bins = np.arange(0., factor, dr)
        rplot = bins+(dr*0.5)

    
    if i == 0.:
        for j in list(range(1, len(ssfr_err) -1 )):
            if (ssfr_err[j] > 10.* ssfr_err[j-1]) & (ssfr_err[j] > 1.):
                ssfr_err[j] = (ssfr_err[j-1] + ssfr_err[j+1]) * 0.5
    elif i == 1:
        for j in list(range(1, len(ssfr_err) -1 )):
            if (np.isnan(ssfr_err[j])) or (ssfr_err[j] > 1.):
                ssfr_err[j] = (ssfr_err[j-1] + ssfr_err[j+1]) * 0.5
    
    ax[0].plot(rplot, ssfr_ratio, color=colors[i], linestyle='-')
    ax[0].fill_between(rplot, ssfr_ratio - ssfr_err, ssfr_ratio + ssfr_err, color=colors[i], alpha=0.25)

ax[0].set_ylim(-2, 2)
ax[0].set_xlim(0, 3.5)
ax[0].set_xlabel(xlabel)
ax[0].set_ylabel(r'$ (\textrm{sSFR}_{\textrm{sat}} - \textrm{sSFR}_{\textrm{cen}} ) / \textrm{sSFR}_{\textrm{cen}}$')

# plot H1:
for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+'s50j7k_star_forming_rotated_faceon_sat_cent_ratio_data.h5') as f:
        h1_ratio = f['h1_ratio_'+b].value
        h1_err = np.abs(f['h1_err_'+b].value)

    ax[1].plot(rplot, h1_ratio, color=colors[i], linestyle='-')
    ax[1].fill_between(rplot, h1_ratio - h1_err, h1_ratio + h1_err, color=colors[i], alpha=0.25)

ax[1].set_ylim(-2, 2)
ax[1].set_xlim(0, 3.5)
ax[1].set_xlabel(xlabel)
ax[1].set_ylabel(r'$ (\Sigma_{\textrm{HI, sat}} - \Sigma_{\textrm{HI, cen}} ) / \Sigma_{\textrm{HI, cen}}$')



plt.savefig(plot_dir+'sat_cent_ratio.png')
plt.clf()
