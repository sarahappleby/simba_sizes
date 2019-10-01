import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from plotting_methods import plot_belfiore
import sys

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 14})

cmap = cm.get_cmap('viridis')
sf_colors = [cmap(0.8), cmap(0.6), cmap(0.3)]
cmap = cm.get_cmap('plasma')
gv_colors = [cmap(0.8), cmap(0.6), cmap(0.25)]

mass_b18 = [10.0,  10.5, 11., 11.5]

data_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/data/'
plot_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/plots/'
bin_labels = ['10.0-10.5', '10.5-11.0', '>11.0']
plot_labels = [r'$10.0 < \textrm{log} (M_* / M_{\odot}) < 10.5$', 
                r'$10.5 < \textrm{log} (M_* / M_{\odot}) < 11.0$', 
                r'$ \textrm{log} (M_* / M_{\odot}) > 11.0$']
xlabel = r'$ R / R_{half}$'

bin_labels.reverse()
gv_colors.reverse()
sf_colors.reverse()
plot_labels.reverse()

fig, ax = plt.subplots(1, 2, figsize=(13, 6))
ax = ax.flatten()

# plot h2:
for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+'sf_rot_h2-weighted_z_data.h5') as f:
        no_gals = f['all_no_gals_'+b].value
        tukey = f['all_tukey_'+b].value
        err = f['all_err_'+b].value
    if i == 0:
        n = len(np.array(tukey))
        dr = 0.2
        factor = dr*n
        bins = np.arange(0., factor, dr)
        rplot = bins+(dr*0.5)

    ax[0].plot(rplot, tukey, color=sf_colors[i], linestyle='-')
    ax[0].fill_between(rplot, tukey - err, tukey + err, color=sf_colors[i], alpha=0.25)

ax[0].set_ylim(-2.2, -1.3)
ax[0].set_xlim(0, 5)
ax[0].set_xlabel(xlabel)
ax[0].set_ylabel(r'$\textrm{log} Z_{H2}$')

for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+'gv_rot_h2-weighted_z_data.h5') as f:
        no_gals = f['all_no_gals_'+b].value
        tukey = f['all_tukey_'+b].value
        err = f['all_err_'+b].value
    if i == 0:
        n = len(np.array(tukey))
        dr = 0.2
        factor = dr*n
        bins = np.arange(0., factor, dr)
        rplot = bins+(dr*0.5)

    ax[1].plot(rplot, tukey, color=sf_colors[i], linestyle='-')
    ax[1].fill_between(rplot, tukey - err, tukey + err, color=sf_colors[i], alpha=0.25)

ax[1].set_ylim(-2.2, -1.3)
ax[1].set_xlim(0, 5)
ax[1].set_xlabel(xlabel)
ax[1].set_ylabel(r'$\textrm{log} Z_{H2}$')


plt.savefig(plot_dir+'h2-weighted_z.png')
plt.show()

