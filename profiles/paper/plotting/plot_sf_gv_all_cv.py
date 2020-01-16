import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from plotting_methods import plot_belfiore, plot_res_limit
import sys

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 14})

softening = 0.5  * 2.8# kpc/h

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

fig, ax = plt.subplots(3, 2, figsize=(13, 15))
ax = ax.flatten()

# plot SFR:
for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+'sf_rot_sfr_data.h5') as f:
        no_gals = f['all_no_gals_'+b].value
        mean_tukey = f['all_jk_'+b].value
        tukey = f['all_tukey_'+b].value
        err = f['all_err_'+b].value
    if i == 0:
        n = len(np.array(tukey))
        dr = 0.2
        factor = dr*n
        bins = np.arange(0., factor, dr)
        rplot = bins+(dr*0.5)
    ax[0].plot(rplot, tukey, color=sf_colors[i], linestyle='-', label=r'$\textrm{SF};\ $' + plot_labels[i])
    ax[0].fill_between(rplot, tukey - err, tukey + err, color=sf_colors[i], alpha=0.25)
 
    with h5py.File(data_dir+'gv_rot_sfr_data.h5') as f:
        no_gals = f['all_no_gals_'+b].value
        mean_tukey = f['all_jk_'+b].value
        tukey = f['all_tukey_'+b].value
        err = f['all_err_'+b].value
    ax[0].plot(rplot, tukey, color=gv_colors[i], linestyle='-')
    ax[0].fill_between(rplot, tukey - err, tukey + err, color=gv_colors[i], alpha=0.25)
ax[0].set_ylim(-5, -1)
ax[0].set_xlim(0, 5)
"""
if gals == 'all':
    ax[0].set_xlim(0, 3.5)
    ax[0].set_ylim(-4.5, -1)
else:
    ax[0].set_xlim(0., 2)
    ax[0].set_ylim(-5.0, -1.5)
"""
ax[0].set_xlabel(xlabel)
ax[0].set_ylabel(r'$ \textrm{log} (\Sigma_{\textrm{SFR}} / M_{\odot}\textrm{yr}^{-1} \textrm{kpc}^{-2})$')
ax[0].legend(fontsize=12)
plot_res_limit(ax[0], softening, 'star_forming', sf_colors)

# plot sSFR:
for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+'sf_rot_ssfr_data.h5') as f:
        no_gals = f['all_no_gals_'+b].value
        mean_tukey = f['all_jk_'+b].value
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
    
    with h5py.File(data_dir+'gv_rot_ssfr_data.h5') as f:
        no_gals = f['all_no_gals_'+b].value
        mean_tukey = f['all_jk_'+b].value
        tukey = f['all_tukey_'+b].value
        err = f['all_err_'+b].value
    ax[1].plot(rplot, tukey, color=gv_colors[i], linestyle='-', label=r'$\textrm{GV};\ $' + plot_labels[i])
    ax[1].fill_between(rplot, tukey - err, tukey + err, color=gv_colors[i], alpha=0.25)
ax[1].set_ylim(-12.5, -9.5)
ax[1].set_xlim(0, 5)
ax[1].set_xlabel(xlabel)
ax[1].set_ylabel(r'$\textrm{log} (\textrm{sSFR} / \textrm{yr}^{-1})$')
ax[1].legend(fontsize=12)
plot_res_limit(ax[1], softening, 'star_forming', sf_colors)

# plot HI:
for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+'sf_rot_h1_data.h5') as f:
        no_gals = f['all_no_gals_'+b].value
        mean_tukey = f['all_jk_'+b].value
        tukey = f['all_tukey_'+b].value
        err = f['all_err_'+b].value
    if i == 0:
        n = len(np.array(tukey))
        dr = 0.2
        factor = dr*n
        bins = np.arange(0., factor, dr)
        rplot = bins+(dr*0.5)
    tukey /= 6.
    err /=6.
    ax[2].plot(rplot, tukey, color=sf_colors[i], linestyle='-')
    ax[2].fill_between(rplot, tukey - err, tukey + err, color=sf_colors[i], alpha=0.25)

    with h5py.File(data_dir+'gv_rot_h1_data.h5') as f:
        no_gals = f['all_no_gals_'+b].value
        mean_tukey = f['all_jk_'+b].value
        tukey = f['all_tukey_'+b].value
        err = f['all_err_'+b].value
    tukey /= 6.
    err /=6.
    ax[2].plot(rplot, tukey, color=gv_colors[i], linestyle='-')
    ax[2].fill_between(rplot, tukey - err, tukey + err, color=gv_colors[i], alpha=0.25)

ax[2].set_ylim(0.7,1.4)
ax[2].set_xlim(0, 5)
ax[2].set_xlabel(xlabel)
ax[2].set_ylabel(r'$ \textrm{log} (\Sigma_{HI} / M_{\odot}\textrm{pc}^{-2})$')
plot_res_limit(ax[2], softening, 'star_forming', sf_colors)

# plot H2:
for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+'sf_rot_h2_data.h5') as f:
        no_gals = f['all_no_gals_'+b].value
        mean_tukey = f['all_jk_'+b].value
        tukey = f['all_tukey_'+b].value
        err = f['all_err_'+b].value
    if i == 0:
        n = len(np.array(tukey))
        dr = 0.2
        factor = dr*n
        bins = np.arange(0., factor, dr)
        rplot = bins+(dr*0.5)
    tukey /= 6.0
    err /=6.
    ax[3].plot(rplot, tukey, color=sf_colors[i], linestyle='-')
    ax[3].fill_between(rplot, tukey - err, tukey + err, color=sf_colors[i], alpha=0.25)

    with h5py.File(data_dir+'gv_rot_h2_data.h5') as f:
        no_gals = f['all_no_gals_'+b].value
        mean_tukey = f['all_jk_'+b].value
        tukey = f['all_tukey_'+b].value
        err = f['all_err_'+b].value
    tukey /= 6.0
    err /=6.
    ax[3].plot(rplot, tukey, color=gv_colors[i], linestyle='-')
    ax[3].fill_between(rplot, tukey - err, tukey + err, color=gv_colors[i], alpha=0.25)

ax[3].axhline(1.1, c='k', ls='--', lw=1.)
ax[3].set_ylim(0.7, 1.4)
ax[3].set_xlim(0, 5)
ax[3].set_xlabel(xlabel)
ax[3].set_ylabel(r'$\textrm{log} (\Sigma_{H_2} / M_{\odot}\textrm{pc}^{-2})$')
plot_res_limit(ax[3], softening, 'star_forming', sf_colors)

# plot sfe:
for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+'sf_rot_sfe_data.h5') as f:
        no_gals = f['all_no_gals_'+b].value
        mean_tukey = f['all_jk_'+b].value
        tukey = f['all_tukey_'+b].value
        err = f['all_err_'+b].value
    if i == 0:
        n = len(np.array(tukey))
        dr = 0.2
        factor = dr*n
        bins = np.arange(0., factor, dr)
        rplot = bins+(dr*0.5)

    ax[4].plot(rplot, tukey, color=sf_colors[i], linestyle='-')
    ax[4].fill_between(rplot, tukey - err, tukey + err, color=sf_colors[i], alpha=0.25)

    with h5py.File(data_dir+'gv_rot_sfe_data.h5') as f:
        no_gals = f['all_no_gals_'+b].value
        mean_tukey = f['all_jk_'+b].value
        tukey = f['all_tukey_'+b].value
        err = f['all_err_'+b].value
    ax[4].plot(rplot, tukey, color=gv_colors[i], linestyle='-')
    ax[4].fill_between(rplot, tukey - err, tukey + err, color=gv_colors[i], alpha=0.25)
ax[4].set_ylim(-11, -9.2)
ax[4].set_xlim(0, 5)
ax[4].set_xlabel(xlabel)
ax[4].set_ylabel(r'$ \textrm{log} (\textrm{SFE} / \textrm{yr}^{-1})$')
plot_res_limit(ax[4], softening, 'star_forming', sf_colors)


# plot fh2:
for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+'sf_rot_fh2_data.h5') as f:
        no_gals = f['all_no_gals_'+b].value
        mean_tukey = f['all_jk_'+b].value
        tukey = f['all_tukey_'+b].value
        err = f['all_err_'+b].value
    if i == 0:
        n = len(np.array(tukey))
        dr = 0.2
        factor = dr*n
        bins = np.arange(0., factor, dr)
        rplot = bins+(dr*0.5)

    ax[5].plot(rplot, tukey, color=sf_colors[i], linestyle='-')
    ax[5].fill_between(rplot, tukey - err, tukey + err, color=sf_colors[i], alpha=0.25)

    with h5py.File(data_dir+'gv_rot_fh2_data.h5') as f:
        no_gals = f['all_no_gals_'+b].value
        mean_tukey = f['all_jk_'+b].value
        tukey = f['all_tukey_'+b].value
        err = f['all_err_'+b].value
    ax[5].plot(rplot, tukey, color=gv_colors[i], linestyle='-')
    ax[5].fill_between(rplot, tukey - err, tukey + err, color=gv_colors[i], alpha=0.25)
ax[5].set_ylim(-3, )
ax[5].set_xlim(0, 5)
ax[5].set_xlabel(xlabel)
ax[5].set_ylabel(r'$\textrm{log} (f_{H2})$')
plot_res_limit(ax[5], softening, 'star_forming', sf_colors)


plt.savefig(plot_dir+'sf_gv_all_cv.png')
plt.show()

