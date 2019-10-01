import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import sys
from plotting_methods import *

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 14})

data_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/data/'
plot_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/plots/'

masks = [2, 3, 4]
bin_labels = ['10.0-10.5', '10.5-11.0', '>11.0']
plot_labels = [r'$10.0 < \textrm{log} (M_* / M_{\odot}) < 10.5$', r'$10.5 < \textrm{log} (M_* / M_{\odot}) < 11.0$', r'$ \textrm{log} (M_* / M_{\odot}) > 11.0$']

cmap = cm.get_cmap('plasma')
colors = [cmap(0.65), cmap(0.5), cmap(0.3)]
colors = [cmap(0.8), cmap(0.6), cmap(0.3)]

cmap = cm.get_cmap('viridis')
colors_sim = [cmap(0.6), cmap(0.35), cmap(0.15)]
colors_sim = [cmap(0.8), cmap(0.6), cmap(0.3)]

wind = sys.argv[1]
angle = 'rot'
mass_b18=[10.0,  10.5, 11., 11.5]

fig, ax = plt.subplots(1, 2, figsize=(15, 6))

# for the star forming galaxies:

for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+wind+'_sf_'+angle+'_ssfr_data.h5') as f:
        no_gals = f['all_no_gals_'+b].value
        tukey = f['all_tukey_'+b].value
        small = f['all_small_scale_'+b].value
        large = f['all_large_scale_'+b].value
        err = f['all_err_'+b].value
    if i == 0:
        n = len(np.array(tukey))
        dr = 0.2
        factor = dr*n
        bins = np.arange(0., factor, dr)
        rplot = bins+(dr*0.5)
    if wind == 's50nojet':
        windlab = r'$\textbf{No-jet}$; '
    elif wind == 's50nox':
        windlab = r'$\textbf{No-Xray}$; '
    ax[0].plot(rplot, tukey, color=colors[i], linestyle='-', label=windlab+ plot_labels[i])
    ax[0].fill_between(rplot, tukey - err, tukey + err, color=colors[i], alpha=0.25)

ax[0].set_xlabel(r'$R / R_{half}$', fontsize=16)
ax[0].set_ylabel(r'$\textrm{log} (\textrm{sSFR} / \textrm{yr}^{-1})$', fontsize=16)
ax[0].set_xlim(0, 2.)
ax[0].set_ylim(-13.5, -9.5)
ax[0].legend(fontsize=12, loc=4)
ax[0].annotate(r'$\textbf{SF}$', xy=(0.15, 0.1), xycoords='axes fraction',size=16)


# for the green valley galaxies:
for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+wind+'_gv_'+angle+'_ssfr_data.h5') as f:
        no_gals = f['all_no_gals_'+b].value
        tukey = f['all_tukey_'+b].value
        small = f['all_small_scale_'+b].value
        large = f['all_large_scale_'+b].value
        err = f['all_cv_err_'+b].value
    if i == 0:
        n = len(np.array(tukey))
        dr = 0.2
        factor = dr*n
        bins = np.arange(0., factor, dr)
        rplot = bins+(dr*0.5)
    isinf = np.where(np.isinf(tukey))[0]
    if len(isinf) > 0.:
        tukey[isinf[0]] = 0.5*(tukey[isinf[0] - 1] + tukey[isinf[0] + 1]) 
        err[isinf[0]] = 0.5*(err[isinf[0] - 1] + err[isinf[0] + 1])

    ax[1].plot(rplot, tukey, color=colors[i], linestyle='-')
    ax[1].fill_between(rplot, tukey - err, tukey + err, color=colors[i], alpha=0.25)

ax[1].set_xlabel(r'$R / R_{half}$', fontsize=16)
ax[1].set_ylabel(r'$\textrm{log} (\textrm{sSFR} / \textrm{yr}^{-1})$', fontsize=16)
ax[1].set_xlim(0, 2.0)
ax[1].set_ylim(-13.5, -9.5)

for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+'sf_'+angle+'_ssfr_data.h5') as f:
        no_gals = f['all_no_gals_'+b].value
        tukey = f['all_tukey_'+b].value
        small = f['all_small_scale_'+b].value
        large = f['all_large_scale_'+b].value
        cv_err = f['all_cv_err_'+b].value
    if i == 0:
        n = len(np.array(tukey))
        dr = 0.2
        factor = dr*n
        bins = np.arange(0., factor, dr)
        rplot = bins+(dr*0.5)
    ax[0].plot(rplot, tukey, color=colors_sim[i], linestyle='-')
    ax[0].fill_between(rplot, tukey - cv_err, tukey + cv_err, color=colors_sim[i], alpha=0.25)

for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+'gv_'+angle+'_ssfr_data.h5') as f:
        no_gals = f['all_no_gals_'+b].value
        tukey = f['all_tukey_'+b].value
        small = f['all_small_scale_'+b].value
        large = f['all_large_scale_'+b].value
        cv_err = f['all_cv_err_'+b].value
    if i == 0:
        n = len(np.array(tukey))
        dr = 0.2
        factor = dr*n
        bins = np.arange(0., factor, dr)
        rplot = bins+(dr*0.5)

    ax[1].plot(rplot, tukey, color=colors_sim[i], linestyle='-', label=r'$\textbf{Simba; }$'+ plot_labels[i])
    ax[1].fill_between(rplot, tukey - cv_err, tukey + cv_err, color=colors_sim[i], alpha=0.25)
ax[1].legend(fontsize=12, loc=1)
ax[1].annotate(r'$\textbf{GV}$', xy=(0.15, 0.9), xycoords='axes fraction',size=16)

plt.savefig(plot_dir+'ssfr_profiles_'+wind+'_belfiore.png')
#plt.show()
#plt.clf()
