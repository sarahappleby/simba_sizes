import h5py
import numpy as np 
import matplotlib.pyplot as plt
import sys

from plotting_methods import tukey_biweight

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 12})

profile_dir = sys.argv[1]
if 'star_forming' in profile_dir:
    mode = 'sf'
elif 'green_valley' in profile_dir:
    mode = 'gv'

results_dir = profile_dir

masks = [2, 3, 4]
bin_labels = [r'$10.0 - 10.5$', r'$10.5 - 11.0$', r'$> 11.0$']

no_gals = np.zeros(3)

for m in masks:

	with h5py.File(profile_dir+'mask_'+str(m)+'_all_profiles.h5', 'r') as f:
    		star_m = f['sm'].value
                gas_m = f['gm'].value
		gas_h1 = f['h1'].value
		gas_h2 = f['h2'].value
        
        if m == 0:
            n = star_m.shape[1]
            dr = 0.2
            factor = dr*n
            bins = np.arange(0., factor, dr)

            h1_frac_tukey = np.zeros((3, n)); h1_frac_large_scale = np.zeros((3, n)); h1_frac_small_scale = np.zeros((3, n))
            h1_mass_tukey = np.zeros((3, n)); h1_mass_large_scale = np.zeros((3, n)); h1_mass_small_scale = np.zeros((3, n))
            h1_mass_sm_median = np.zeros((3, n)); h1_mass_sm_lower = np.zeros((3, n)); h1_mass_sm_higher = np.zeros((3, n)) 
            h1_mass_sm_mean = np.zeros((3, n)); h1_mass_sm_sigma = np.zeros((3, n))

            h2_frac_tukey = np.zeros((3, n)); h2_frac_large_scale = np.zeros((3, n)); h2_frac_small_scale = np.zeros((3, n))
            h2_mass_tukey = np.zeros((3, n)); h2_mass_large_scale = np.zeros((3, n)); h2_mass_small_scale = np.zeros((3, n))
            h2_mass_sm_median = np.zeros((3, n)); h2_mass_sm_lower = np.zeros((3, n)); h2_mass_sm_higher = np.zeros((3, n))
            h2_mass_sm_mean = np.zeros((3, n)); h2_mass_sm_sigma = np.zeros((3, n))


    	h1_mass_sm = gas_h1*gas_m / star_m
        h2_mass_sm = gas_h2*gas_m / star_m

	no_gals[m] = len(star_m)

        
        # hi and h2:

        tukey, scale = tukey_biweight(gas_h1)
	h1_frac_tukey[m] = np.log10(tukey)
        h1_frac_large_scale[m] = scale / (np.log(10.)*tukey)
        h1_frac_small_scale[m] = scale / (np.sqrt(no_gals[m])* np.log(10.)*tukey)
        tukey, scale = tukey_biweight(gas_h1*gas_m)
        h1_mass_tukey[m] = np.log10(tukey)
        h1_mass_large_scale[m] = scale / (np.log(10.)*tukey)
        h1_mass_small_scale[m] = scale / (np.sqrt(no_gals[m])* np.log(10.)*tukey)
        
        h1_mass_sm_median[m] = np.nanpercentile(h1_mass_sm, '50', axis=0)
        h1_mass_sm_lower[m] = np.nanpercentile(h1_mass_sm, '25', axis=0)
        h1_mass_sm_higher[m] = np.nanpercentile(h1_mass_sm, '75', axis=0)
        h1_mass_sm_mean[m] = np.mean(h1_mass_sm, axis=0)
        h1_mass_sm_sigma[m] = np.std(h1_mass_sm, axis=0)


        tukey, scale = tukey_biweight(gas_h2)
        h2_frac_tukey[m] = np.log10(tukey)
        h2_frac_large_scale[m] = scale / (np.log(10.)*tukey)
        h2_frac_small_scale[m] = scale / (np.sqrt(no_gals[m])* np.log(10.)*tukey)
        tukey, scale = tukey_biweight(gas_h2*gas_m)
        h2_mass_tukey[m] = np.log10(tukey)
        h2_mass_large_scale[m] = scale / (np.log(10.)*tukey)
        h2_mass_small_scale[m] = scale / (np.sqrt(no_gals[m])* np.log(10.)*tukey)

        h2_mass_sm_median[m] = np.nanpercentile(h2_mass_sm, '50', axis=0)
        h2_mass_sm_lower[m] = np.nanpercentile(h2_mass_sm, '25', axis=0)
        h2_mass_sm_higher[m] = np.nanpercentile(h2_mass_sm, '75', axis=0)
        h2_mass_sm_mean[m] = np.mean(h2_mass_sm, axis=0)
        h2_mass_sm_sigma[m] = np.std(h2_mass_sm, axis=0)

for m in range(len(bin_labels)):
    plt.plot(bins+(dr*0.5), h1_frac_tukey[m], marker='.', markersize=4, linestyle='-', label=bin_labels[m] +', '+str(int(no_gals[m]))+' galaxies')
    #plt.fill_between(bins+(dr*0.5), h1_frac_tukey[m] - h1_frac_large_scale[m], h1_frac_tukey[m] + h1_frac_large_scale[m], alpha=0.1)
    plt.fill_between(bins+(dr*0.5), h1_frac_tukey[m] - h1_frac_small_scale[m], h1_frac_tukey[m] + h1_frac_small_scale[m], alpha=0.3)
plt.legend()
plt.xlabel(r'$R_{half}$')
plt.ylabel(r'$\Sigma_{f\ HI} (\textrm{kpc}^{-2})$')
if factor == 2:
    plt.xlim(0, 1.5)
else:
    plt.xlim(0, factor)
plt.savefig(results_dir+'h1_frac_tukey.png')
plt.clf()

for m in range(len(bin_labels)):
    plt.plot(bins+(dr*0.5), h1_mass_tukey[m], marker='.', markersize=4, linestyle='-', label=bin_labels[m] +', '+str(int(no_gals[m]))+' galaxies')
    #plt.fill_between(bins+(dr*0.5), h1_mass_tukey[m] - h1_mass_large_scale[m], h1_mass_tukey[m] + h1_mass_large_scale[m], alpha=0.1)
    plt.fill_between(bins+(dr*0.5), h1_mass_tukey[m] - h1_mass_small_scale[m], h1_mass_tukey[m] + h1_mass_small_scale[m], alpha=0.3)
plt.legend()
plt.xlabel(r'$R_{half}$')
plt.ylabel(r'$\Sigma_{HI} (M_{\odot}\textrm{kpc}^{-2})$')
if factor == 2:
    plt.xlim(0, 1.5)
else:
    plt.xlim(0, factor)
plt.savefig(results_dir+'h1_mass_tukey.png')
plt.clf()

for m in range(len(bin_labels)):
        plt.plot(bins+(dr*0.5), h1_mass_sm_median[m], marker='.', markersize=4, linestyle='--', label=bin_labels[m] +', '+str(int(no_gals[m]))+' galaxies')
        plt.fill_between(bins+(dr*0.5), h1_mass_sm_lower[m], h1_mass_sm_higher[m], alpha=0.3)
plt.legend()
plt.xlabel(r'$R_{half}$')
plt.ylabel(r'$ M_{HI} / M_*$')
plt.savefig(results_dir+'h1_mass_sm_medians.png')
plt.clf()
for m in range(len(bin_labels)):
        plt.errorbar(bins+(dr*0.5), h1_mass_sm_mean[m], yerr=h1_mass_sm_sigma[m], marker='.', markersize=4, linestyle='--',
                                label=bin_labels[m] +', '+str(int(no_gals[m]))+' galaxies')
plt.legend()
plt.xlabel('$R_{half}$')
plt.ylabel(r'$ M_{HI} / M_*$')
plt.savefig(results_dir+'h1_mass_sm_means.png')
plt.clf()


for m in range(len(bin_labels)):
    plt.plot(bins+(dr*0.5), h2_frac_tukey[m], marker='.', markersize=4, linestyle='-', label=bin_labels[m] +', '+str(int(no_gals[m]))+' galaxies')
    #plt.fill_between(bins+(dr*0.5), h2_frac_tukey[m] - h2_frac_large_scale[m], h2_frac_tukey[m] + h2_frac_large_scale[m], alpha=0.1)
    plt.fill_between(bins+(dr*0.5), h2_frac_tukey[m] - h2_frac_small_scale[m], h2_frac_tukey[m] + h2_frac_small_scale[m], alpha=0.3)
plt.legend()
plt.xlabel(r'$R_{half}$')
plt.ylabel(r'$\Sigma_{f\ H_2} (\textrm{kpc}^{-2})$')
if factor == 2:
    plt.xlim(0, 1.5)
else:
    plt.xlim(0, factor)
plt.savefig(results_dir+'h2_frac_tukey.png')
plt.clf()

for m in range(len(bin_labels)):
    plt.plot(bins+(dr*0.5), h2_mass_tukey[m], marker='.', markersize=4, linestyle='-', label=bin_labels[m] +', '+str(int(no_gals[m]))+' galaxies')
    #plt.fill_between(bins+(dr*0.5), h2_mass_tukey[m] - h2_mass_large_scale[m], h2_mass_tukey[m] + h2_mass_large_scale[m], alpha=0.1)
    plt.fill_between(bins+(dr*0.5), h2_mass_tukey[m] - h2_mass_small_scale[m], h2_mass_tukey[m] + h2_mass_small_scale[m], alpha=0.3)
plt.legend()
plt.xlabel(r'$R_{half}$')
plt.ylabel(r'$\Sigma_{H_2} (M_{\odot}\textrm{kpc}^{-2})$')
if factor == 2:
    plt.xlim(0, 1.5)
else:
    plt.xlim(0, factor)
plt.savefig(results_dir+'h2_mass_tukey.png')
plt.clf()

for m in range(len(bin_labels)):
        plt.plot(bins+(dr*0.5), h2_mass_sm_median[m], marker='.', markersize=4, linestyle='--', label=bin_labels[m] +', '+str(int(no_gals[m]))+' galaxies')
        plt.fill_between(bins+(dr*0.5), h2_mass_sm_lower[m], h2_mass_sm_higher[m], alpha=0.3)
plt.legend()
plt.xlabel(r'$R_{half}$')
plt.ylabel(r'$ M_{H_2} / M_*$')
plt.savefig(results_dir+'h2_mass_sm_medians.png')
plt.clf()
for m in range(len(bin_labels)):
        plt.errorbar(bins+(dr*0.5), h2_mass_sm_mean[m], yerr=h2_mass_sm_sigma[m], marker='.', markersize=4, linestyle='--',
                                label=bin_labels[m] +', '+str(int(no_gals[m]))+' galaxies')
plt.legend()
plt.xlabel('$R_{half}$')
plt.ylabel(r'$ M_{H_2} / M_*$')
plt.savefig(results_dir+'h2_mass_sm_means.png')
plt.clf()

