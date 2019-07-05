import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys
from plotting_methods import *

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 14})

model = sys.argv[1]
wind = sys.argv[2]
snap = sys.argv[3]

basic_dir = '/home/sapple/simba_sizes/profiles/ssfr/extended_profiles/'
gv_centrals_dir = basic_dir + 'centrals/'+model+'_'+snap+'/'+wind+'/green_valley/random_orientation/'
gv_sats_dir = basic_dir + 'satellites/'+model+'_'+snap+'/'+wind+'/green_valley/random_orientation/'
sf_centrals_dir = basic_dir + 'centrals/'+model+'_'+snap+'/'+wind+'/star_forming/random_orientation/'
sf_sats_dir = basic_dir + 'satellites/'+model+'_'+snap+'/'+wind+'/star_forming/random_orientation/'
results_dir = '/home/sapple/simba_sizes/profiles/plotting/plots/'

bin_labels_centrals = [r'$10.0 < \textrm{log} (M_* / M_{\odot}) < 10.5$', r'$10.5 < \textrm{log} (M_* / M_{\odot}) < 11.0$', r'$ \textrm{log} (M_* / M_{\odot}) > 11.0$']
bin_labels_sats = [r'$10.0 < \textrm{log} (M_* / M_{\odot}) < 10.5$', r'$\textrm{log} (M_* / M_{\odot}) > 10.5$']
masks = [2, 3, 4]
colors = ['b', 'm', 'r']
mass_b18=[10.0,  10.5, 11.0, 11.5]


fig, ax = plt.subplots(2, 2, figsize=(15, 15))

# CENTRALS:
# for the star forming galaxies:
no_gals = np.zeros(len(masks))
for i, m in enumerate(masks):

    with h5py.File(sf_centrals_dir+'mask_'+str(m)+'_all_profiles.h5', 'r') as f:
        star_m = f['sm'].value
        gas_sfr = f['gas_sfr'].value
    """
    if m == 3:
        with h5py.File(sf_centrals_dir+'mask_'+str(m+1)+'_all_profiles.h5', 'r') as f:
            star_m = np.concatenate((star_m, f['sm'].value))
            gas_sfr = np.concatenate((gas_sfr, f['gas_sfr'].value))
    """

    if i == 0:
        n = star_m.shape[1]
        dr = 0.2
        factor = dr*n
        bins = np.arange(0., factor, dr)
        gas_ssfr_tukey = np.zeros((len(masks), n)); gas_ssfr_large_scale = np.zeros((len(masks), n)); gas_ssfr_small_scale = np.zeros((len(masks), n))

    no_gals[i] = len(gas_sfr)
    gas_ssfr = gas_sfr / star_m
    tukey, scale = tukey_biweight(gas_ssfr)
    gas_ssfr_tukey[i] = np.log10(tukey)
    gas_ssfr_large_scale[i] = scale / (np.log(10.)*tukey)
    gas_ssfr_small_scale[i] = scale / (np.sqrt(no_gals[i])* np.log(10.)*tukey)

plot_belfiore(ax[0][0], 'sf', colors, mass_b18=mass_b18)
for m in range(len(bin_labels_centrals)):
    ax[0][0].plot(bins+(dr*0.5), gas_ssfr_tukey[m], color=colors[m], marker='.', markersize=4, linestyle='--', label=bin_labels_centrals[m] +'; '+str(int(no_gals[m]))+' galaxies')
    if m == 0:
        ax[0][0].fill_between(bins+(dr*0.5), gas_ssfr_tukey[m] - gas_ssfr_large_scale[m], gas_ssfr_tukey[m] + gas_ssfr_large_scale[m], color=colors[m], alpha=0.1)
    ax[0][0].fill_between(bins+(dr*0.5), gas_ssfr_tukey[m] - gas_ssfr_small_scale[m], gas_ssfr_tukey[m] + gas_ssfr_small_scale[m], color=colors[m], alpha=0.3)

ax[0][0].set_xlim(0., 1.5)
ax[0][0].set_ylim(-12.5, -9.5)
ax[0][0].set_xlabel(r'$R_{half}$', fontsize=16)
ax[0][0].set_ylabel(r'$\textrm{log} (\textrm{sSFR} / \textrm{yr}^{-1})$', fontsize=16)
ax[0][0].legend()

# for the green valley galaxies:
no_gals = np.zeros(len(masks))
for i, m in enumerate(masks):

    with h5py.File(gv_centrals_dir+'mask_'+str(m)+'_all_profiles.h5', 'r') as f:
        star_m = f['sm'].value
        gas_sfr = f['gas_sfr'].value

    """
    if m == 3:
        with h5py.File(gv_centrals_dir+'mask_'+str(m+1)+'_all_profiles.h5', 'r') as f:
            star_m = np.concatenate((star_m, f['sm'].value))
            gas_sfr = np.concatenate((gas_sfr, f['gas_sfr'].value))
    """

    if i == 0:
        n = star_m.shape[1]
        dr = 0.2
        factor = dr*n
        bins = np.arange(0., factor, dr)
        gas_ssfr_tukey = np.zeros((len(masks), n)); gas_ssfr_large_scale = np.zeros((len(masks), n)); gas_ssfr_small_scale = np.zeros((len(masks), n))

    no_gals[i] = len(gas_sfr)
    gas_ssfr = gas_sfr / star_m
    tukey, scale = tukey_biweight(gas_ssfr)
    gas_ssfr_tukey[i] = np.log10(tukey)
    gas_ssfr_large_scale[i] = scale / (np.log(10.)*tukey)
    gas_ssfr_small_scale[i] = scale / (np.sqrt(no_gals[i])* np.log(10.)*tukey)

plot_belfiore(ax[0][1], 'gv', colors, mass_b18=mass_b18)
for m in range(len(bin_labels_centrals)):
    ax[0][1].plot(bins+(dr*0.5), gas_ssfr_tukey[m], color=colors[m], marker='.', markersize=4, linestyle='--', label=str(int(no_gals[m]))+' galaxies')
    if m == 0:
        ax[0][1].fill_between(bins+(dr*0.5), gas_ssfr_tukey[m] - gas_ssfr_large_scale[m], gas_ssfr_tukey[m] + gas_ssfr_large_scale[m], color=colors[m], alpha=0.1)
    ax[0][1].fill_between(bins+(dr*0.5), gas_ssfr_tukey[m] - gas_ssfr_small_scale[m], gas_ssfr_tukey[m] + gas_ssfr_small_scale[m], color=colors[m], alpha=0.3)

ax[0][1].set_xlim(0, 1.5)
ax[0][1].set_ylim(-12.5, -9.5)
ax[0][1].set_xlabel(r'$R_{half}$', fontsize=16)
ax[0][1].set_ylabel(r'$\textrm{log} (\textrm{sSFR} / \textrm{yr}^{-1})$', fontsize=16)
ax[0][1].legend()


# SATELLITES:
# For the star forming galaxies:
no_gals = np.zeros(len(masks))
for i, m in enumerate(masks[:-1]):

    with h5py.File(sf_sats_dir+'mask_'+str(m)+'_all_profiles.h5', 'r') as f:
        star_m = f['sm'].value
        gas_sfr = f['gas_sfr'].value

    if m == 3:
        with h5py.File(sf_sats_dir+'mask_'+str(m+1)+'_all_profiles.h5', 'r') as f:
            star_m = np.concatenate((star_m, f['sm'].value))
            gas_sfr = np.concatenate((gas_sfr, f['gas_sfr'].value))

    if i == 0:
        n = star_m.shape[1]
        dr = 0.2
        factor = dr*n
        bins = np.arange(0., factor, dr)
        gas_ssfr_tukey = np.zeros((len(masks), n)); gas_ssfr_large_scale = np.zeros((len(masks), n)); gas_ssfr_small_scale = np.zeros((len(masks), n))

    no_gals[i] = len(gas_sfr)
    gas_ssfr = gas_sfr / star_m
    tukey, scale = tukey_biweight(gas_ssfr)
    gas_ssfr_tukey[i] = np.log10(tukey)
    gas_ssfr_large_scale[i] = scale / (np.log(10.)*tukey)
    gas_ssfr_small_scale[i] = scale / (np.sqrt(no_gals[i])* np.log(10.)*tukey)

for m in range(len(bin_labels_sats)):
    ax[1][0].plot(bins+(dr*0.5), gas_ssfr_tukey[m], color=colors[m], marker='.', markersize=4, linestyle=':', label=bin_labels_sats[m]+'; '+str(int(no_gals[m]))+' galaxies')
    if m == 0:
        ax[1][0].fill_between(bins+(dr*0.5), gas_ssfr_tukey[m] - gas_ssfr_large_scale[m], gas_ssfr_tukey[m] + gas_ssfr_large_scale[m], color=colors[m], alpha=0.1)
    ax[1][0].fill_between(bins+(dr*0.5), gas_ssfr_tukey[m] - gas_ssfr_small_scale[m], gas_ssfr_tukey[m] + gas_ssfr_small_scale[m], color=colors[m], alpha=0.3)

ax[1][0].set_xlim(0, 1.5)
ax[1][0].set_ylim(-12.5, -9.5)
ax[1][0].set_xlabel(r'$R_{half}$', fontsize=16)
ax[1][0].set_ylabel(r'$\textrm{log} (\textrm{sSFR} / \textrm{yr}^{-1})$', fontsize=16)
ax[1][0].legend()

# for the green valley galaxies:
no_gals = np.zeros(len(masks))
for i, m in enumerate(masks[:-1]):

    with h5py.File(gv_sats_dir+'mask_'+str(m)+'_all_profiles.h5', 'r') as f:
        star_m = f['sm'].value
        gas_sfr = f['gas_sfr'].value

    if m == 3:
        with h5py.File(gv_sats_dir+'mask_'+str(m+1)+'_all_profiles.h5', 'r') as f:
            star_m = np.concatenate((star_m, f['sm'].value))
            gas_sfr = np.concatenate((gas_sfr, f['gas_sfr'].value))

    if i == 0:
        n = star_m.shape[1]
        dr = 0.2
        factor = dr*n
        bins = np.arange(0., factor, dr)
        gas_ssfr_tukey = np.zeros((len(masks), n)); gas_ssfr_large_scale = np.zeros((len(masks), n)); gas_ssfr_small_scale = np.zeros((len(masks), n))

    no_gals[i] = len(gas_sfr)
    gas_ssfr = gas_sfr / star_m
    tukey, scale = tukey_biweight(gas_ssfr)
    gas_ssfr_tukey[i] = np.log10(tukey)
    gas_ssfr_large_scale[i] = scale / (np.log(10.)*tukey)
    gas_ssfr_small_scale[i] = scale / (np.sqrt(no_gals[i])* np.log(10.)*tukey)

for m in range(len(bin_labels_sats)):
    ax[1][1].plot(bins+(dr*0.5), gas_ssfr_tukey[m], color=colors[m], marker='.', markersize=4, linestyle=':', label=str(int(no_gals[m]))+' galaxies')
    if m == 0:
        ax[1][1].fill_between(bins+(dr*0.5), gas_ssfr_tukey[m] - gas_ssfr_large_scale[m], gas_ssfr_tukey[m] + gas_ssfr_large_scale[m], color=colors[m], alpha=0.1)
    ax[1][1].fill_between(bins+(dr*0.5), gas_ssfr_tukey[m] - gas_ssfr_small_scale[m], gas_ssfr_tukey[m] + gas_ssfr_small_scale[m], color=colors[m], alpha=0.3)

ax[1][1].set_xlabel(r'$R_{half}$', fontsize=16)
ax[1][1].set_ylabel(r'$\textrm{log} (\textrm{sSFR} / \textrm{yr}^{-1})$', fontsize=16)
ax[1][1].set_xlim(0, 1.5)
ax[1][1].set_ylim(-12.5, -9.5)
ax[1][1].legend()

plt.savefig(results_dir+'ssfr_centrals_satellites_'+wind+'.png')
plt.clf()

