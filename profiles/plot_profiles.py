import h5py
import numpy as np 
import matplotlib.pyplot as plt
import sys

profile_dir = sys.argv[1]

results_dir = profile_dir

mass_bins = [10., 10.5, 11.]
bin_labels = ['10.0 - 10.5', '10.5 - 11.0', '> 11.0']

n = 10
factor = 2.
dr = factor / n
bins = np.arange(0., factor, dr)

sm_median = np.zeros((3, n)); sm_lower = np.zeros((3, n)); sm_higher = np.zeros((3, n)); sm_mean = np.zeros((3, n)); sm_sigma = np.zeros((3, n))
star_sfr_median = np.zeros((3, n)); star_sfr_lower = np.zeros((3, n)); star_sfr_higher = np.zeros((3, n)); star_sfr_mean = np.zeros((3, n)); star_sfr_sigma = np.zeros((3, n))
star_ssfr_median = np.zeros((3, n)); star_ssfr_lower = np.zeros((3, n)); star_ssfr_higher = np.zeros((3, n)); star_ssfr_mean = np.zeros((3, n)); star_ssfr_sigma = np.zeros((3, n))

gas_sfr_median = np.zeros((3, n)); gas_sfr_lower = np.zeros((3, n)); gas_sfr_higher = np.zeros((3, n)); gas_sfr_mean = np.zeros((3, n)); gas_sfr_sigma = np.zeros((3, n))
gas_ssfr_median = np.zeros((3, n)); gas_ssfr_lower = np.zeros((3, n)); gas_ssfr_higher = np.zeros((3, n)); gas_ssfr_mean = np.zeros((3, n)); gas_ssfr_sigma = np.zeros((3, n))
gas_h1_median = np.zeros((3, n)); gas_h1_lower = np.zeros((3, n)); gas_h1_higher = np.zeros((3, n)); gas_h1_mean = np.zeros((3, n)); gas_h1_sigma = np.zeros((3, n))
gas_h2_median = np.zeros((3, n)); gas_h2_lower = np.zeros((3, n)); gas_h2_higher = np.zeros((3, n)); gas_h2_mean = np.zeros((3, n)); gas_h2_sigma = np.zeros((3, n))

no_gals = np.zeros(3)

for m in range(3):

	with h5py.File(profile_dir+'mask_'+str(m)+'_all_profiles.h5', 'r') as f:
		use_star_sfr = f['star_sfr'].value
		use_star_m = f['sm'].value
		use_gas_sfr = f['gas_sfr'].value
		use_gas_h1 = f['h1'].value
		use_gas_h2 = f['h2'].value

	no_gals[m] = len(use_star_m)

	sm_median[m] = np.percentile(use_star_m, '50', axis=0)
	sm_lower[m] = np.percentile(use_star_m, '25', axis=0)
	sm_higher[m] = np.percentile(use_star_m, '75', axis=0)
	sm_mean[m] = np.mean(use_star_m, axis=0)
	sm_sigma[m] = np.std(use_star_m, axis=0)

	star_sfr_median[m] = np.nanpercentile(use_star_sfr, '50', axis=0)
	star_sfr_lower[m] = np.nanpercentile(use_star_sfr, '25', axis=0)
	star_sfr_higher[m] = np.nanpercentile(use_star_sfr, '75', axis=0)
	star_sfr_mean[m] = np.mean(use_star_sfr, axis=0)
	star_sfr_sigma[m] = np.std(use_star_sfr, axis=0)

	# in manga paper, ssfr profiles are ratio of median sfr and median sm

	star_ssfr_median[m] = star_sfr_median[m] / sm_median[m]
	star_ssfr_lower[m] = star_sfr_lower[m] / sm_lower[m]
	star_ssfr_higher[m] = star_sfr_higher[m] / sm_higher[m]
	star_ssfr_mean[m] = star_sfr_mean[m] / sm_mean[m]
	star_ssfr_sigma[m] = star_ssfr_mean[m] * ((star_sfr_sigma[m] / star_sfr_mean[m])**2 + (sm_sigma[m] / sm_mean[m])**2)**0.5

	star_ssfr_median[m][star_ssfr_median[m] == 0.] = 1.e-20
	star_ssfr_lower[m][star_ssfr_lower[m] == 0.] = 1.e-20
	star_ssfr_higher[m][star_ssfr_higher[m] == 0.] = 1.e-20
	star_ssfr_mean[m][star_ssfr_mean[m] == 0.] = 1.e-20
	star_ssfr_sigma[m] /= (star_ssfr_mean[m]*np.log(10.))


	gas_sfr_median[m] = np.nanpercentile(use_gas_sfr, '50', axis=0)
	gas_sfr_lower[m] = np.nanpercentile(use_gas_sfr, '25', axis=0)
	gas_sfr_higher[m] = np.nanpercentile(use_gas_sfr, '75', axis=0)
	gas_sfr_mean[m] = np.mean(use_gas_sfr, axis=0)
	gas_sfr_sigma[m] = np.std(use_gas_sfr, axis=0)

	gas_ssfr_median[m] = gas_sfr_median[m] / sm_median[m]
	gas_ssfr_lower[m] = gas_sfr_lower[m] / sm_lower[m]
	gas_ssfr_higher[m] = gas_sfr_higher[m] / sm_higher[m]
	gas_ssfr_mean[m] = gas_sfr_mean[m] / sm_mean[m]
	gas_ssfr_sigma[m] = gas_ssfr_mean[m] * ((gas_sfr_sigma[m] / gas_sfr_mean[m])**2 + (sm_sigma[m] / sm_mean[m])**2)**0.5

	gas_ssfr_median[m][gas_ssfr_median[m] == 0.] = 1.e-20
	gas_ssfr_lower[m][gas_ssfr_lower[m] == 0.] = 1.e-20
	gas_ssfr_higher[m][gas_ssfr_higher[m] == 0.] = 1.e-20
	gas_ssfr_mean[m][gas_ssfr_mean[m] == 0.] = 1.e-20
	gas_ssfr_sigma[m] /= (gas_ssfr_mean[m]*np.log(10.))

	gas_h1_median[m] = np.nanpercentile(use_gas_h1, '50', axis=0)
	gas_h1_lower[m] = np.nanpercentile(use_gas_h1, '25', axis=0)
	gas_h1_higher[m] = np.nanpercentile(use_gas_h1, '75', axis=0)
	gas_h1_mean[m] = np.nanmean(use_gas_h1, axis=0)
	gas_h1_sigma[m] = np.nanstd(use_gas_h1, axis=0)

	gas_h2_median[m] = np.nanpercentile(use_gas_h2, '50', axis=0)
	gas_h2_lower[m] = np.nanpercentile(use_gas_h2, '25', axis=0)
	gas_h2_higher[m] = np.nanpercentile(use_gas_h2, '75', axis=0)
	gas_h2_mean[m] = np.nanmean(use_gas_h2, axis=0)
	gas_h2_sigma[m] = np.nanstd(use_gas_h2, axis=0)

for m in range(len(mass_bins)):
	plt.semilogy(bins+(dr*0.5), sm_median[m], marker='.', markersize=4, linestyle='--', label=bin_labels[m] +', '+str(no_gals[m])+' galaxies')
	plt.fill_between(bins+(dr*0.5), sm_lower[m], sm_higher[m], alpha=0.3)
plt.legend()
plt.xlabel('R half *')
plt.ylabel('M* surface density (Msun/R half bin^2)')
plt.xlim(0, 2)
plt.savefig(results_dir+'sm_medians.png')
plt.clf()
for m in range(len(mass_bins)):
	plt.errorbar(bins+(dr*0.5), sm_mean[m], yerr=sm_sigma[m], marker='.', markersize=4, linestyle='--', 
				label=bin_labels[m] +', '+str(no_gals[m])+' galaxies')
plt.yscale('log')
plt.legend()
plt.xlabel('R half *')
plt.ylabel('M* surface density (Msun/R half bin^2)')
plt.xlim(0, 2)
plt.savefig(results_dir+'sm_means.png')
plt.clf()


for m in range(len(mass_bins)):
	plt.plot(bins+(dr*0.5), star_sfr_median[m], marker='.', markersize=4, linestyle='--', label=bin_labels[m] +', '+str(no_gals[m])+' galaxies')
	plt.fill_between(bins+(dr*0.5), star_sfr_lower[m], star_sfr_higher[m], alpha=0.3)
plt.legend()
plt.xlabel('R half *')
plt.ylabel('SFR surface density (Msun/yr /R half bin^2)')
plt.xlim(0, 2)
plt.savefig(results_dir+'star_sfr_medians.png')
plt.clf()
for m in range(len(mass_bins)):
	plt.errorbar(bins+(dr*0.5), star_sfr_mean[m], yerr=star_sfr_sigma[m], marker='.', markersize=4, linestyle='--', 
				label=bin_labels[m] +', '+str(no_gals[m])+' galaxies')
plt.legend()
plt.xlabel('R half *')
plt.ylabel('SFR surface density (Msun/yr /R half bin^2)')
plt.xlim(0, 2)
plt.savefig(results_dir+'star_sfr_means.png')
plt.clf()


for m in range(len(mass_bins)):
	plt.plot(bins+(dr*0.5), np.log10(star_ssfr_median[m]), marker='.', markersize=4, linestyle='--', label=bin_labels[m] +', '+str(no_gals[m])+' galaxies')
	plt.fill_between(bins+(dr*0.5), np.log10(star_ssfr_lower[m]), np.log10(star_ssfr_higher[m]), alpha=0.3)
plt.legend()
plt.xlabel('R half *')
plt.ylabel('log sSFR (yr^-1)')
plt.xlim(0, 2)
plt.ylim(-13, )
plt.savefig(results_dir+'star_ssfr_medians.png')
plt.clf()
for m in range(len(mass_bins)):
	plt.errorbar(bins+(dr*0.5), np.log10(star_ssfr_mean[m]), yerr=star_ssfr_sigma[m], marker='.', markersize=4, linestyle='--', 
				label=bin_labels[m] +', '+str(no_gals[m])+' galaxies')
plt.legend()
plt.xlabel('R half *')
plt.ylabel('log sSFR (yr^-1)')
plt.xlim(0, 2)
plt.ylim(-13, )
plt.savefig(results_dir+'star_ssfr_means.png')
plt.clf()


for m in range(len(mass_bins)):
	plt.plot(bins+(dr*0.5), gas_sfr_median[m], marker='.', markersize=4, linestyle='--', label=bin_labels[m] +', '+str(no_gals[m])+' galaxies')
	plt.fill_between(bins+(dr*0.5), gas_sfr_lower[m], gas_sfr_higher[m], alpha=0.3)
plt.legend()
plt.xlabel('R half *')
plt.ylabel('SFR surface density (Msun/yr /R half bin^2)')
plt.xlim(0, 2)
plt.ylim(0, )
plt.savefig(results_dir+'gas_sfr_medians.png')
plt.clf()
for m in range(len(mass_bins)):
	plt.errorbar(bins+(dr*0.5), gas_sfr_mean[m], yerr=gas_sfr_sigma[m], marker='.', markersize=4, linestyle='--', 
				label=bin_labels[m] +', '+str(no_gals[m])+' galaxies')
plt.legend()
plt.xlabel('R half *')
plt.ylabel('SFR surface density (Msun/yr /R half bin^2)')
plt.xlim(0, 2)
plt.ylim(0, )
plt.savefig(results_dir+'gas_sfr_means.png')
plt.clf()


for m in range(len(mass_bins)):
	plt.plot(bins+(dr*0.5), np.log10(gas_ssfr_median[m]), marker='.', markersize=4, linestyle='--', label=bin_labels[m] +', '+str(no_gals[m])+' galaxies')
	plt.fill_between(bins+(dr*0.5), np.log10(gas_ssfr_lower[m]), np.log10(gas_ssfr_higher[m]), alpha=0.3)
plt.legend()
plt.xlabel('R half *')
plt.ylabel('log sSFR (yr^-1)')
plt.xlim(0, 2)
plt.ylim(-13, )
plt.savefig(results_dir+'gas_ssfr_medians.png')
plt.clf()
for m in range(len(mass_bins)):
	plt.errorbar(bins+(dr*0.5), np.log10(gas_ssfr_mean[m]), yerr=gas_ssfr_sigma[m], marker='.', markersize=4, linestyle='--', 
				label=bin_labels[m] +', '+str(no_gals[m])+' galaxies')
plt.legend()
plt.xlabel('R half *')
plt.ylabel('log sSFR (yr^-1)')
plt.xlim(0, 2)
plt.ylim(-13, )
plt.savefig(results_dir+'gas_ssfr_means.png')
plt.clf()


for m in range(len(mass_bins)):
	plt.plot(bins+(dr*0.5), gas_h1_median[m], marker='.', markersize=4, linestyle='--', label=bin_labels[m] +', '+str(no_gals[m])+' galaxies')
	plt.fill_between(bins+(dr*0.5), gas_h1_lower[m], gas_h1_higher[m], alpha=0.3)
plt.legend()
plt.xlabel('R half *')
plt.ylabel('HI surface density (Msun /R half bin^2)')
plt.xlim(0, 2)
plt.ylim(0, )
plt.savefig(results_dir+'h1_medians.png')
plt.clf()
for m in range(len(mass_bins)):
	plt.errorbar(bins+(dr*0.5), gas_h1_mean[m], yerr=gas_h1_sigma[m], marker='.', markersize=4, linestyle='--', 
				label=bin_labels[m] +', '+str(no_gals[m])+' galaxies')
plt.legend()
plt.xlabel('R half *')
plt.ylabel('HI surface density (Msun /R half bin^2)')
plt.xlim(0, 2)
plt.ylim(0, )
plt.savefig(results_dir+'h1_means.png')
plt.clf()


for m in range(len(mass_bins)):
	plt.plot(bins+(dr*0.5), gas_h2_median[m], marker='.', markersize=4, linestyle='--', label=bin_labels[m] +', '+str(no_gals[m])+' galaxies')
	plt.fill_between(bins+(dr*0.5), gas_h2_lower[m], gas_h2_higher[m], alpha=0.3)
plt.legend()
plt.xlabel('R half *')
plt.ylabel('HII surface density (Msun /R half bin^2)')
plt.xlim(0, 2)
plt.ylim(0, )
plt.savefig(results_dir+'h2_medians.png')
plt.clf()
for m in range(len(mass_bins)):
	plt.errorbar(bins+(dr*0.5), gas_h2_mean[m], yerr=gas_h2_sigma[m], marker='.', markersize=4, linestyle='--', 
				label=bin_labels[m] +', '+str(no_gals[m])+' galaxies')
plt.legend()
plt.xlabel('R half *')
plt.ylabel('HII surface density (Msun /R half bin^2)')
plt.xlim(0, 2)
plt.ylim(0, )
plt.savefig(results_dir+'h2_means.png')
plt.clf()