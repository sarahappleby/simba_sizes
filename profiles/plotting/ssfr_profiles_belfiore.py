import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys
from plotting_methods import *

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 12})

model = sys.argv[1]
wind = sys.argv[2]
snap = sys.argv[3]

gv_dir = '/home/sapple/simba_sizes/profiles/ssfr/extended_profiles/centrals/'+model+'_'+snap+'/'+wind+'/green_valley/random_orientation/'
sf_dir = '/home/sapple/simba_sizes/profiles/ssfr/extended_profiles/centrals/'+model+'_'+snap+'/'+wind+'/star_forming/random_orientation/'

#old:
#gv_dir = '/home/sapple/simba_sizes/profiles/ssfr/inner_galaxy/m100n1024_145/green_valley/random_orientation/bh_centered/'
#sf_dir = '/home/sapple/simba_sizes/profiles/ssfr/inner_galaxy/m100n1024_145/star_forming/random_orientation/bh_centered/'
results_dir = '/home/sapple/simba_sizes/profiles/plotting/'

masks = [2, 3, 4]
bin_labels = [r'$10.0 - 10.5$', r'$10.5 - 11.0$', r'$> 11.0$']
colors = ['b', 'm', 'r']
mass_c=[10.0,  10.5, 11., 11.5]

fig, ax = plt.subplots(1, 2, figsize=(15, 6))

no_gals = np.zeros(3)
# for the star forming galaxies:
for i, m in enumerate(masks):

    with h5py.File(sf_dir+'mask_'+str(m)+'_all_profiles.h5', 'r') as f:
        star_m = f['sm'].value
        gas_sfr = f['gas_sfr'].value

    if i == 0:
        n = star_m.shape[1]
        dr = 0.2
        factor = dr*n
        bins = np.arange(0., factor, dr)
        gas_ssfr_tukey = np.zeros((3, n)); gas_ssfr_large_scale = np.zeros((3, n)); gas_ssfr_small_scale = np.zeros((3, n))
    
    no_gals[i] = len(gas_sfr)
    gas_ssfr = gas_sfr / star_m
    tukey, scale = tukey_biweight(gas_ssfr)
    gas_ssfr_tukey[i] = np.log10(tukey)
    gas_ssfr_large_scale[i] = scale / (np.log(10.)*tukey)
    gas_ssfr_small_scale[i] = scale / (np.sqrt(no_gals[i])* np.log(10.)*tukey)

plot_belfiore(ax[0], 'sf', colors, mass_c=mass_c)
for m in range(len(bin_labels)):
    ax[0].plot(bins+(dr*0.5), gas_ssfr_tukey[m], color=colors[m], marker='.', markersize=4, linestyle='--', label=bin_labels[m] +', '+str(int(no_gals[m]))+' galaxies')
    if m == 0:
        ax[0].fill_between(bins+(dr*0.5), gas_ssfr_tukey[m] - gas_ssfr_large_scale[m], gas_ssfr_tukey[m] + gas_ssfr_large_scale[m], color=colors[m], alpha=0.1)
    ax[0].fill_between(bins+(dr*0.5), gas_ssfr_tukey[m] - gas_ssfr_small_scale[m], gas_ssfr_tukey[m] + gas_ssfr_small_scale[m], color=colors[m], alpha=0.3)
ax[0].set_xlabel(r'$R_{half}$', fontsize=16)
ax[0].set_ylabel(r'$\textrm{log} (\textrm{sSFR} / \textrm{yr}^{-1})$', fontsize=16)
ax[0].set_xlim(0, 1.5)
ax[0].set_ylim(-12.5, -9.)
ax[0].legend(fontsize=10)

no_gals = np.zeros(3)
# for the green valley galaxies:
for i, m in enumerate(masks):

    with h5py.File(gv_dir+'mask_'+str(m)+'_all_profiles.h5', 'r') as f:
        star_m = f['sm'].value
        gas_sfr = f['gas_sfr'].value
       
    if i == 0:
        n = star_m.shape[1]
        dr = 0.2 
        factor = dr*n
        bins = np.arange(0., factor, dr) 
        gas_ssfr_tukey = np.zeros((3, n)); gas_ssfr_large_scale = np.zeros((3, n)); gas_ssfr_small_scale = np.zeros((3, n))
        
    no_gals[i] = len(gas_sfr)
    gas_ssfr = gas_sfr / star_m
    tukey, scale = tukey_biweight(gas_ssfr)
    gas_ssfr_tukey[i] = np.log10(tukey)
    gas_ssfr_large_scale[i] = scale / (np.log(10.)*tukey)
    gas_ssfr_small_scale[i] = scale / (np.sqrt(no_gals[i])* np.log(10.)*tukey)

plot_belfiore(ax[1], 'gv', colors, mass_c=mass_c)
for m in range(len(bin_labels)):
    ax[1].plot(bins+(dr*0.5), gas_ssfr_tukey[m], color=colors[m], marker='.', markersize=4, linestyle='--', label=bin_labels[m] +', '+str(int(no_gals[m]))+' galaxies')
    if m == 0:
        ax[1].fill_between(bins+(dr*0.5), gas_ssfr_tukey[m] - gas_ssfr_large_scale[m], gas_ssfr_tukey[m] + gas_ssfr_large_scale[m], color=colors[m], alpha=0.1)
    ax[1].fill_between(bins+(dr*0.5), gas_ssfr_tukey[m] - gas_ssfr_small_scale[m], gas_ssfr_tukey[m] + gas_ssfr_small_scale[m], color=colors[m], alpha=0.3)
ax[1].set_xlabel(r'$R_{half}$', fontsize=16)
ax[1].set_ylabel(r'$\textrm{log} (\textrm{sSFR} / \textrm{yr}^{-1})$', fontsize=16)
ax[1].set_xlim(0, 1.5)
ax[1].set_ylim(-12.5, -9.)
ax[1].legend(fontsize=10)

plt.savefig(results_dir+'ssfr_profiles_'+wind+'_belfiore.png')
plt.clf()
