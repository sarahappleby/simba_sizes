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
gals = 'all'
angle = 'rot'
mass_b18 = [10.0,  10.5, 11., 11.5]

data_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/data/'
plot_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/plots/'
bin_labels = ['10.0-10.5', '10.5-11.0', '>11.0']
plot_labels = [r'$10.0 < \textrm{log} (M_* / M_{\odot}) < 10.5$', 
                r'$10.5 < \textrm{log} (M_* / M_{\odot}) < 11.0$', 
                r'$ \textrm{log} (M_* / M_{\odot}) > 11.0$']
xlabel = r'$ R / R_{half}$'

fig, ax = plt.subplots(1, 2, figsize=(13, 6))
ax = ax.flatten()

# plot SFE:
for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+selection+'_'+angle+'_sfe_data.h5') as f:
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

if selection == 'sf':
    ax[0].set_ylim(-11, -9)
    ax[0].set_xlim(0, 3)
else:
    ax[0].set_ylim(-11, -9)
    ax[0].set_xlim(0, 1.5)
ax[0].set_xlabel(xlabel)
ax[0].set_ylabel(r'$ \textrm{log} (\textrm{SFE} / M_{\odot}\textrm{yr}^{-1})$')

# plot fh2:
for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+selection+'_'+angle+'_fh2_data.h5') as f:
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

    ax[1].plot(rplot, tukey, color=colors[i], linestyle='-')
    if i == 1:
        ax[1].fill_between(rplot, tukey - large, tukey + large, color=colors[i], alpha=0.1)
    ax[1].fill_between(rplot, tukey - small, tukey + small, color=colors[i], alpha=0.3)

ax[1].set_ylim(-3, )
ax[1].set_xlabel(xlabel)
ax[1].set_ylabel(r'$\textrm{log} (f_{H2})$')

plt.savefig(plot_dir+selection+'_rot_sfe_fh2_radial.png')
plt.show()
plt.clf()
