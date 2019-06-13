import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys
from plotting_methods import *

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 12})

wind = 's50j7k'

basic_dir = '/home/sapple/simba_sizes/profiles/ssfr/extended_profiles/'
gv_centrals_dir = basic_dir + 'centrals/m100n1024_151/'+wind+'/green_valley/random_orientation/'
gv_sats_dir = basic_dir + 'satellites/m100n1024_151/'+wind+'/green_valley/random_orientation/'
sf_centrals_dir = basic_dir + 'centrals/m100n1024_151/'+wind+'/star_forming/random_orientation/'
sf_sats_dir = basic_dir + 'satellites/m100n1024_151/'+wind+'/star_forming/random_orientation/'

bin_labels = [r'$9.0 - 9.5$', r'$9.5 - 10.0$', r'$10.0 - 10.5$', r'$>10.5$']
masks = [0, 1, 2, 3]

colors = ['o', 'g', 'b', 'm', 'r']

fig, ax = plt.subplots(1, 2, figsize=(15, 6))

# for the star forming galaxies:
# CENTRALS:
no_gals = np.zeros(len(masks))
for i, m in enumerate(masks):

    with h5py.File(sf_centrals_dir+'mask_'+str(m)+'_all_profiles.h5', 'r') as f:
        star_m = f['sm'].value
        gas_sfr = f['gas_sfr'].value

    if m == 3:
        with h5py.File(sf_centrals_dir+'mask_'+str(m+1)+'_all_profiles.h5', 'r') as f:
            star_m = np.concatenate((star_m, f['sm'].value))
            gas_m = np.concatenate((gas_m, f['gm'].value))
            gas_h1 = np.concatenate((gas_h1, f['h1'].value))
            gas_h2 = np.concatenate((gas_h2, f['h2'].value))

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

for m in range(len(bin_labels)):
    ax[0].plot(bins+(dr*0.5), gas_ssfr_tukey[m], color=colors[m], marker='.', markersize=4, linestyle='--', label=bin_labels[m] +', '+str(int(no_gals[m]))+' galaxies')
    if m == 0:
        ax[0].fill_between(bins+(dr*0.5), gas_ssfr_tukey[m] - gas_ssfr_large_scale[m], gas_ssfr_tukey[m] + gas_ssfr_large_scale[m], color=colors[m], alpha=0.1)
    ax[0].fill_between(bins+(dr*0.5), gas_ssfr_tukey[m] - gas_ssfr_small_scale[m], gas_ssfr_tukey[m] + gas_ssfr_small_scale[m], color=colors[m], alpha=0.3)

# SATELLITES:
no_gals = np.zeros(len(masks))
for i, m in enumerate(masks):

    with h5py.File(sf_sats_dir+'mask_'+str(m)+'_all_profiles.h5', 'r') as f:
        star_m = f['sm'].value
        gas_sfr = f['gas_sfr'].value

    if m == 3:
        with h5py.File(sf_sats_dir+'mask_'+str(m+1)+'_all_profiles.h5', 'r') as f:
            star_m = np.concatenate((star_m, f['sm'].value))
            gas_m = np.concatenate((gas_m, f['gm'].value))
            gas_h1 = np.concatenate((gas_h1, f['h1'].value))
            gas_h2 = np.concatenate((gas_h2, f['h2'].value))

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

for m in range(len(bin_labels)):
    ax[0].plot(bins+(dr*0.5), gas_ssfr_tukey[m], color=colors[m], marker='.', markersize=4, linestyle=':', label=bin_labels[m] +', '+str(int(no_gals[m]))+' galaxies')
    if m == 0:
        ax[0].fill_between(bins+(dr*0.5), gas_ssfr_tukey[m] - gas_ssfr_large_scale[m], gas_ssfr_tukey[m] + gas_ssfr_large_scale[m], color=colors[m], alpha=0.1)
    ax[0].fill_between(bins+(dr*0.5), gas_ssfr_tukey[m] - gas_ssfr_small_scale[m], gas_ssfr_tukey[m] + gas_ssfr_small_scale[m], color=colors[m], alpha=0.3)

ax[0].set_xlabel(r'$R_{half}$')
ax[0].set_ylabel(r'$\textrm{log} (\textrm{sSFR} / \textrm{yr}^{-1})$')
ax[0].set_xlim(0, 1.5)
ax[0].set_ylim(-12.5, -9.)
ax[0].legend()

# for the green valley galaxies:
#CENTRALS:
no_gals = np.zeros(len(masks))
for i, m in enumerate(masks):

    with h5py.File(gf_centrals_dir+'mask_'+str(m)+'_all_profiles.h5', 'r') as f:
        star_m = f['sm'].value
        gas_sfr = f['gas_sfr'].value

    if m == 3:
        with h5py.File(gf_centrals_dir+'mask_'+str(m+1)+'_all_profiles.h5', 'r') as f:
            star_m = np.concatenate((star_m, f['sm'].value))
            gas_m = np.concatenate((gas_m, f['gm'].value))
            gas_h1 = np.concatenate((gas_h1, f['h1'].value))
            gas_h2 = np.concatenate((gas_h2, f['h2'].value))

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

for m in range(len(bin_labels)):
    ax[1].plot(bins+(dr*0.5), gas_ssfr_tukey[m], color=colors[m], marker='.', markersize=4, linestyle='--', label=bin_labels[m] +', '+str(int(no_gals[m]))+' galaxies')
    if m == 0:
        ax[1].fill_between(bins+(dr*0.5), gas_ssfr_tukey[m] - gas_ssfr_large_scale[m], gas_ssfr_tukey[m] + gas_ssfr_large_scale[m], color=colors[m], alpha=0.1)
    ax[1].fill_between(bins+(dr*0.5), gas_ssfr_tukey[m] - gas_ssfr_small_scale[m], gas_ssfr_tukey[m] + gas_ssfr_small_scale[m], color=colors[m], alpha=0.3)

#SATELLITES
no_gals = np.zeros(len(masks))
for i, m in enumerate(masks):

    with h5py.File(gf_sats_dir+'mask_'+str(m)+'_all_profiles.h5', 'r') as f:
        star_m = f['sm'].value
        gas_sfr = f['gas_sfr'].value

    if m == 3:
        with h5py.File(gf_sats_dir+'mask_'+str(m+1)+'_all_profiles.h5', 'r') as f:
            star_m = np.concatenate((star_m, f['sm'].value))
            gas_m = np.concatenate((gas_m, f['gm'].value))
            gas_h1 = np.concatenate((gas_h1, f['h1'].value))
            gas_h2 = np.concatenate((gas_h2, f['h2'].value))

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

for m in range(len(bin_labels)):
    ax[1].plot(bins+(dr*0.5), gas_ssfr_tukey[m], color=colors[m], marker='.', markersize=4, linestyle=':', label=bin_labels[m] +', '+str(int(no_gals[m]))+' galaxies')
    if m == 0:
        ax[1].fill_between(bins+(dr*0.5), gas_ssfr_tukey[m] - gas_ssfr_large_scale[m], gas_ssfr_tukey[m] + gas_ssfr_large_scale[m], color=colors[m], alpha=0.1)
    ax[1].fill_between(bins+(dr*0.5), gas_ssfr_tukey[m] - gas_ssfr_small_scale[m], gas_ssfr_tukey[m] + gas_ssfr_small_scale[m], color=colors[m], alpha=0.3)

ax[1].set_xlabel(r'$R_{half}$')
ax[1].set_ylabel(r'$\textrm{log} (\textrm{sSFR} / \textrm{yr}^{-1})$')
ax[1].set_xlim(0, 1.5)
ax[1].set_ylim(-12.5, -9.)
ax[1].legend()

plt.savefig(results_dir+'ssfr_centrals_satellites.png')
plt.clf()

