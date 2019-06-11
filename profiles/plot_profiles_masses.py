import h5py
import numpy as np 
import matplotlib.pyplot as plt
import sys

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 12})

def tukey_biweight(x, c=5.0, epsilon=1.e-11):
    median = np.nanpercentile(x, '50', axis=0)
    mad = np.nanmedian(np.abs(x - median), axis=0)
    u = (x - median)/(c*mad + epsilon)
    weights = np.zeros(u.shape)
    mask = np.abs(u) < 1.
    weights[mask] = (1. - u[mask]**2)**2
    tukey = np.nansum(x*weights, axis=0) / np.sum(weights, axis=0)
    n = x.shape[0]
    num = np.sqrt(np.nansum(((x - tukey)**2. * (1. - u**2.)**4.)*mask, axis=0))
    den = np.abs(np.nansum(((1. - u )*(1 - 5.*(u**2.)))*mask, axis=0))
    scale = np.sqrt(n)*num / den
    t = 1.96 # 97.5% of t distribution with df = max(0.7(n-1), 1)
    return tukey, scale


profile_dir = sys.argv[1]
if 'star_forming' in profile_dir:
    mode = 'sf'
elif 'green_valley' in profile_dir:
    mode = 'gv'

results_dir = profile_dir

#bin_labels = ['10.0 - 10.5', '10.5 - 11.0']

masks = [2, 3, 4]
bin_labels = [r'$10.0 - 10.5$', r'$10.5 - 11.0$', r'$> 11.0$']

no_gals = np.zeros(3)

for m in masks:

	with h5py.File(profile_dir+'mask_'+str(m)+'_all_profiles.h5', 'r') as f:
    		star_m = f['sm'].value
                gas_m = f['gm'].value
                gas_sfr = f['gas_sfr'].value
		gas_h1 = f['h1'].value
		gas_h2 = f['h2'].value
        
        if m == 0:
            n = star_m.shape[1]
            dr = 0.2
            factor = dr*n
            bins = np.arange(0., factor, dr)

            sm_median = np.zeros((3, n)); sm_lower = np.zeros((3, n)); sm_higher = np.zeros((3, n)); sm_mean = np.zeros((3, n)); sm_sigma = np.zeros((3, n))
            
            gas_sfr_median = np.zeros((3, n)); gas_sfr_lower = np.zeros((3, n)); gas_sfr_higher = np.zeros((3, n)); gas_sfr_mean = np.zeros((3, n)); gas_sfr_sigma = np.zeros((3, n))
            gas_ssfr_tukey = np.zeros((3, n)); gas_ssfr_large_scale = np.zeros((3, n)); gas_ssfr_small_scale = np.zeros((3, n))
            gas_ssfr_median = np.zeros((3, n)); gas_ssfr_lower = np.zeros((3, n)); gas_ssfr_higher = np.zeros((3, n)); gas_ssfr_mean = np.zeros((3, n)); gas_ssfr_sigma = np.zeros((3, n))
            
            h1_frac_tukey = np.zeros((3, n)); h1_frac_large_scale = np.zeros((3, n)); h1_frac_small_scale = np.zeros((3, n))
            h1_mass_tukey = np.zeros((3, n)); h1_mass_large_scale = np.zeros((3, n)); h1_mass_small_scale = np.zeros((3, n))
            
            h2_frac_tukey = np.zeros((3, n)); h2_frac_large_scale = np.zeros((3, n)); h2_frac_small_scale = np.zeros((3, n))
            h2_mass_tukey = np.zeros((3, n)); h2_mass_large_scale = np.zeros((3, n)); h2_mass_small_scale = np.zeros((3, n))

    	gas_ssfr = gas_sfr / star_m

	no_gals[m] = len(star_m)

        # star particle profiles

	sm_median[m] = np.percentile(star_m, '50', axis=0)
	sm_lower[m] = np.percentile(star_m, '25', axis=0)
	sm_higher[m] = np.percentile(star_m, '75', axis=0)
	sm_mean[m] = np.mean(star_m, axis=0)
	sm_sigma[m] = np.std(star_m, axis=0)


        # gas particle profiles

	gas_sfr_median[m] = np.nanpercentile(gas_sfr, '50', axis=0)
	gas_sfr_lower[m] = np.nanpercentile(gas_sfr, '25', axis=0)
	gas_sfr_higher[m] = np.nanpercentile(gas_sfr, '75', axis=0)
	gas_sfr_mean[m] = np.mean(gas_sfr, axis=0)
	gas_sfr_sigma[m] = np.std(gas_sfr, axis=0)
       
        tukey, scale = tukey_biweight(gas_ssfr)
        gas_ssfr_tukey[m] = np.log10(tukey)
        gas_ssfr_large_scale[m] = scale / (np.log(10.)*tukey)
        gas_ssfr_small_scale[m] = scale / (np.sqrt(no_gals[m])* np.log(10.)*tukey)
        

        use_gas_ssfr = gas_ssfr.copy()
        use_gas_ssfr[use_gas_ssfr == 0.] = 1.e-14
        use_gas_ssfr = np.log10(use_gas_ssfr)
        gas_ssfr_median[m] = np.nanpercentile(use_gas_ssfr, '50', axis=0)
	gas_ssfr_lower[m] = np.nanpercentile(use_gas_ssfr, '25', axis=0)
	gas_ssfr_higher[m] = np.nanpercentile(use_gas_ssfr, '75', axis=0)
	gas_ssfr_mean[m] = np.mean(use_gas_ssfr, axis=0)
	gas_ssfr_sigma[m] = np.std(use_gas_ssfr, axis=0)
	
        # hi and h2:

        tukey, scale = tukey_biweight(gas_h1)
	h1_frac_tukey[m] = np.log10(tukey)
        h1_frac_large_scale[m] = scale / (np.log(10.)*tukey)
        h1_frac_small_scale[m] = scale / (np.sqrt(no_gals[m])* np.log(10.)*tukey)
        tukey, scale = tukey_biweight(gas_h1*gas_m)
        h1_mass_tukey[m] = np.log10(tukey)
        h1_mass_large_scale[m] = scale / (np.log(10.)*tukey)
        h1_mass_small_scale[m] = scale / (np.sqrt(no_gals[m])* np.log(10.)*tukey)

        tukey, scale = tukey_biweight(gas_h2)
        h2_frac_tukey[m] = np.log10(tukey)
        h2_frac_large_scale[m] = scale / (np.log(10.)*tukey)
        h2_frac_small_scale[m] = scale / (np.sqrt(no_gals[m])* np.log(10.)*tukey)
        tukey, scale = tukey_biweight(gas_h2*gas_m)
        h2_mass_tukey[m] = np.log10(tukey)
        h2_mass_large_scale[m] = scale / (np.log(10.)*tukey)
        h2_mass_small_scale[m] = scale / (np.sqrt(no_gals[m])* np.log(10.)*tukey)

for m in range(len(bin_labels)):
	plt.semilogy(bins+(dr*0.5), sm_median[m], marker='.', markersize=4, linestyle='--', label=bin_labels[m] +', '+str(int(no_gals[m]))+' galaxies')
	plt.fill_between(bins+(dr*0.5), sm_lower[m], sm_higher[m], alpha=0.3)
plt.legend()
plt.xlabel(r'$R_{half}$')
plt.ylabel(r'$\Sigma_{M*} (M_{\odot}\textrm{kpc}^{-2})$')
plt.xlim(0, factor)
plt.savefig(results_dir+'sm_medians.png')
plt.clf()
for m in range(len(bin_labels)):
	plt.errorbar(bins+(dr*0.5), sm_mean[m], yerr=sm_sigma[m], marker='.', markersize=4, linestyle='--', 
				label=bin_labels[m] +', '+str(int(no_gals[m]))+' galaxies')
plt.yscale('log')
plt.legend()
plt.xlabel(r'$R_{half}$')
plt.ylabel(r'$\Sigma_{M*} (M_{\odot}\textrm{kpc}^{-2})$')
plt.xlim(0, factor)
plt.savefig(results_dir+'sm_means.png')
plt.clf()

for m in range(len(bin_labels)):
	plt.plot(bins+(dr*0.5), gas_sfr_median[m], marker='.', markersize=4, linestyle='--', label=bin_labels[m] +', '+str(int(no_gals[m]))+' galaxies')
	plt.fill_between(bins+(dr*0.5), gas_sfr_lower[m], gas_sfr_higher[m], alpha=0.3)
plt.legend()
plt.xlabel(r'$R_{half}$')
plt.ylabel(r'$\Sigma_{\textrm{SFR}} (M_{\odot}\textrm{yr}^{-1} \textrm{kpc}^{-2})$')
plt.xlim(0, factor)
plt.ylim(0, )
plt.savefig(results_dir+'gas_sfr_medians.png')
plt.clf()
for m in range(len(bin_labels)):
	plt.errorbar(bins+(dr*0.5), gas_sfr_mean[m], yerr=gas_sfr_sigma[m], marker='.', markersize=4, linestyle='--', 
				label=bin_labels[m] +', '+str(int(no_gals[m]))+' galaxies')
plt.legend()
plt.xlabel(r'$R_{half}$')
plt.ylabel(r'$\Sigma_{\textrm{SFR}} (M_{\odot}\textrm{yr}^{-1} \textrm{kpc}^{-2})$')
plt.xlim(0, factor)
plt.ylim(0, )
plt.savefig(results_dir+'gas_sfr_means.png')
plt.clf()

for m in range(len(bin_labels)):
    plt.plot(bins+(dr*0.5), gas_ssfr_tukey[m], marker='.', markersize=4, linestyle='-', label=bin_labels[m] +', '+str(int(no_gals[m]))+' galaxies')
    #plt.fill_between(bins+(dr*0.5), gas_ssfr_tukey[m] - gas_ssfr_large_scale[m], gas_ssfr_tukey[m] + gas_ssfr_large_scale[m], alpha=0.1)
    plt.fill_between(bins+(dr*0.5), gas_ssfr_tukey[m] - gas_ssfr_small_scale[m], gas_ssfr_tukey[m] + gas_ssfr_small_scale[m], alpha=0.3)
if mode == 'sf':
    plt.legend(loc=4)
elif mode == 'gv':
    plt.legend(loc=1)
plt.xlabel(r'$R_{half}$')
plt.ylabel(r'$\textrm{log} (\textrm{sSFR} / \textrm{yr}^{-1})$')
if factor == 2:
    plt.xlim(0, 1.5)
else:
    plt.xlim(0, factor)
plt.ylim(-12.5, -9.5)
plt.savefig(results_dir+'gas_ssfr_tukey.png')
plt.clf()
for m in range(len(bin_labels)):
	plt.plot(bins+(dr*0.5), gas_ssfr_median[m], marker='.', markersize=4, linestyle='--', label=bin_labels[m] +', '+str(int(no_gals[m]))+' galaxies')
	plt.fill_between(bins+(dr*0.5), gas_ssfr_lower[m], gas_ssfr_higher[m], alpha=0.3)
plt.legend()
plt.xlabel(r'$R_{half}$')
plt.ylabel(r'$\textrm{log} (\textrm{sSFR} / \textrm{yr}^{-1})$')
if factor == 2:
    plt.xlim(0, 1.5)
else:
    plt.xlim(0, factor)
plt.ylim(-12.5, )
plt.savefig(results_dir+'gas_ssfr_medians.png')
plt.clf()
for m in range(len(bin_labels)):
	plt.errorbar(bins+(dr*0.5), gas_ssfr_mean[m], yerr=gas_ssfr_sigma[m], marker='.', markersize=4, linestyle='--', 
				label=bin_labels[m] +', '+str(int(no_gals[m]))+' galaxies')
plt.legend()
plt.xlabel('$R_{half}$')
plt.ylabel(r'$\textrm{log} (\textrm{sSFR} / \textrm{yr}^{-1})$')
if factor == 2:
    plt.xlim(0, 1.5)
plt.ylim(-12.5, )
plt.savefig(results_dir+'gas_ssfr_means.png')
plt.clf()

for m in range(len(bin_labels)):
    plt.plot(bins+(dr*0.5), h1_frac_tukey[m], marker='.', markersize=4, linestyle='-', label=bin_labels[m] +', '+str(int(no_gals[m]))+' galaxies')
    #plt.fill_between(bins+(dr*0.5), h1_frac_tukey[m] - h1_frac_large_scale[m], gas_ssfr_tukey[m] + gas_ssfr_large_scale[m], alpha=0.1)
    plt.fill_between(bins+(dr*0.5), h1_frac_tukey[m] - h1_frac_small_scale[m], gas_ssfr_tukey[m] + gas_ssfr_small_scale[m], alpha=0.3)
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
    #plt.fill_between(bins+(dr*0.5), h1_mass_tukey[m] - h1_mass_large_scale[m], gas_ssfr_tukey[m] + gas_ssfr_large_scale[m], alpha=0.1)
    plt.fill_between(bins+(dr*0.5), h1_mass_tukey[m] - h1_mass_small_scale[m], gas_ssfr_tukey[m] + gas_ssfr_small_scale[m], alpha=0.3)
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
    plt.plot(bins+(dr*0.5), h2_frac_tukey[m], marker='.', markersize=4, linestyle='-', label=bin_labels[m] +', '+str(int(no_gals[m]))+' galaxies')
    #plt.fill_between(bins+(dr*0.5), h2_frac_tukey[m] - h2_frac_large_scale[m], gas_ssfr_tukey[m] + gas_ssfr_large_scale[m], alpha=0.1)
    plt.fill_between(bins+(dr*0.5), h2_frac_tukey[m] - h2_frac_small_scale[m], gas_ssfr_tukey[m] + gas_ssfr_small_scale[m], alpha=0.3)
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
    #plt.fill_between(bins+(dr*0.5), h2_mass_tukey[m] - h2_mass_large_scale[m], gas_ssfr_tukey[m] + gas_ssfr_large_scale[m], alpha=0.1)
    plt.fill_between(bins+(dr*0.5), h2_mass_tukey[m] - h2_mass_small_scale[m], gas_ssfr_tukey[m] + gas_ssfr_small_scale[m], alpha=0.3)
plt.legend()
plt.xlabel(r'$R_{half}$')
plt.ylabel(r'$\Sigma_{H_2} (M_{\odot}\textrm{kpc}^{-2})$')
if factor == 2:
    plt.xlim(0, 1.5)
else:
    plt.xlim(0, factor)
plt.savefig(results_dir+'h2_mass_tukey.png')
plt.clf()

