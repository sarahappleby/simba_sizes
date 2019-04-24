import h5py
import numpy as np 
import matplotlib.pyplot as plt
import sys

selection = 'green_valley'
bh_center = True
rotate_galaxies = False

data_dir = '/home/sapple/simba_sizes/profiles/science/'
results_dir = data_dir + 'wind_comparison/'

winds = ['s50j7k', 's50nojet', 's50nox']
wind_labels = ['Simba', 'No jet', 'No Xray']

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

for w, wind in enumerate(winds):

	profile_dir = data_dir + wind + '/selection_1/m50n512_145/' + selection + '/'
	if rotate_galaxies:
		results_dir += '/rotated_faceon'
	else:
		results_dir += '/random_orientation'
	if bh_center:
		results_dir += '/bh_centered/'
	else:
		results_dir += '/scom_centered/'

	for m in range(3):
		with h5py.File(profile_dir+'mask_'+str(m)+'_all_profiles.h5', 'r') as f:
		if m == 0:
				use_star_sfr = f['star_sfr'].value
				use_star_m = f['sm'].value
				use_gas_sfr = f['gas_sfr'].value
				use_gas_h1 = f['h1'].value
				use_gas_h2 = f['h2'].value
		else:

			use_star_sfr = np.concatenate((use_star_sfr, f['star_sfr'].value))
			use_star_m = np.concatenate((use_star_m, f['sm'].value))
			use_gas_sfr = np.concatenate((use_gas_sfr, f['gas_sfr'].value))
			use_gas_h1 = np.concatenate((use_gas_h1, f['h1'].value))
			use_gas_h2 = np.concatenate((use_gas_h2, f['h2'].value))

	no_gals[w] = len(use_star_m)


	sm_median[w] = np.percentile(use_star_m, '50', axis=0)
	sm_lower[w] = np.percentile(use_star_m, '25', axis=0)
	sm_higher[w] = np.percentile(use_star_m, '75', axis=0)
	sm_mean[w] = np.mean(use_star_m, axis=0)
	sm_sigma[w] = np.std(use_star_m, axis=0)

	star_sfr_median[w] = np.nanpercentile(use_star_sfr, '50', axis=0)
	star_sfr_lower[w] = np.nanpercentile(use_star_sfr, '25', axis=0)
	star_sfr_higher[w] = np.nanpercentile(use_star_sfr, '75', axis=0)
	star_sfr_mean[w] = np.mean(use_star_sfr, axis=0)
	star_sfr_sigma[w] = np.std(use_star_sfr, axis=0)

	star_ssfr_median[w] = star_sfr_median[w] / sm_median[w]
	star_ssfr_lower[w] = star_sfr_lower[w] / sm_lower[w]
	star_ssfr_higher[w] = star_sfr_higher[w] / sm_higher[w]
	star_ssfr_mean[w] = star_sfr_mean[w] / sm_mean[w]
	star_ssfr_sigma[w] = star_ssfr_mean[w] * ((star_sfr_sigma[w] / star_sfr_mean[w])**2 + (sm_sigma[w] / sm_mean[w])**2)**0.5

	star_ssfr_median[w][star_ssfr_median[w] == 0.] = 1.e-20
	star_ssfr_lower[w][star_ssfr_lower[w] == 0.] = 1.e-20
	star_ssfr_higher[w][star_ssfr_higher[w] == 0.] = 1.e-20
	star_ssfr_mean[w][star_ssfr_mean[w] == 0.] = 1.e-20
	star_ssfr_sigma[w] /= (star_ssfr_mean[w]*np.log(10.))


	gas_sfr_median[w] = np.nanpercentile(use_gas_sfr, '50', axis=0)
	gas_sfr_lower[w] = np.nanpercentile(use_gas_sfr, '25', axis=0)
	gas_sfr_higher[w] = np.nanpercentile(use_gas_sfr, '75', axis=0)
	gas_sfr_mean[w] = np.mean(use_gas_sfr, axis=0)
	gas_sfr_sigma[w] = np.std(use_gas_sfr, axis=0)

	gas_ssfr_median[w] = gas_sfr_median[w] / sm_median[w]
	gas_ssfr_lower[w] = gas_sfr_lower[w] / sm_lower[w]
	gas_ssfr_higher[w] = gas_sfr_higher[w] / sm_higher[w]
	gas_ssfr_mean[w] = gas_sfr_mean[w] / sm_mean[w]
	gas_ssfr_sigma[w] = gas_ssfr_mean[w] * ((gas_sfr_sigma[w] / gas_sfr_mean[w])**2 + (sm_sigma[w] / sm_mean[w])**2)**0.5

	gas_ssfr_median[w][gas_ssfr_median[w] == 0.] = 1.e-20
	gas_ssfr_lower[w][gas_ssfr_lower[w] == 0.] = 1.e-20
	gas_ssfr_higher[w][gas_ssfr_higher[w] == 0.] = 1.e-20
	gas_ssfr_mean[w][gas_ssfr_mean[w] == 0.] = 1.e-20
	gas_ssfr_sigma[w] /= (gas_ssfr_mean[w]*np.log(10.))

	gas_h1_median[w] = np.nanpercentile(use_gas_h1, '50', axis=0)
	gas_h1_lower[w] = np.nanpercentile(use_gas_h1, '25', axis=0)
	gas_h1_higher[w] = np.nanpercentile(use_gas_h1, '75', axis=0)
	gas_h1_mean[w] = np.nanmean(use_gas_h1, axis=0)
	gas_h1_sigma[w] = np.nanstd(use_gas_h1, axis=0)

	gas_h2_median[w] = np.nanpercentile(use_gas_h2, '50', axis=0)
	gas_h2_lower[w] = np.nanpercentile(use_gas_h2, '25', axis=0)
	gas_h2_higher[w] = np.nanpercentile(use_gas_h2, '75', axis=0)
	gas_h2_mean[w] = np.nanmean(use_gas_h2, axis=0)
	gas_h2_sigma[w] = np.nanstd(use_gas_h2, axis=0)



for m in range(len(winds)):
	plt.semilogy(bins+(dr*0.5), sm_median[w], marker='.', markersize=4, linestyle='--', label=wind_labels[w] +', '+str(no_gals[w])+' galaxies')
	plt.fill_between(bins+(dr*0.5), sm_lower[w], sm_higher[w], alpha=0.3)
plt.legend()
plt.xlabel('R half *')
plt.ylabel('M* surface density (Msun/kpc^2)')
plt.xlim(0, 2)
plt.savefig(results_dir+'sm_medians.png')
plt.clf()
for m in range(len(winds)):
	plt.errorbar(bins+(dr*0.5), sm_mean[w], yerr=sm_sigma[m], marker='.', markersize=4, linestyle='--', 
				label=wind_labels[m] +', '+str(no_gals[m])+' galaxies')
plt.yscale('log')
plt.legend()
plt.xlabel('R half *')
plt.ylabel('M* surface density (Msun/kpc^2)')
plt.xlim(0, 2)
plt.savefig(results_dir+'sm_means.png')
plt.clf()


for m in range(len(winds)):
	plt.plot(bins+(dr*0.5), star_sfr_median[m], marker='.', markersize=4, linestyle='--', label=wind_labels[m] +', '+str(no_gals[m])+' galaxies')
	plt.fill_between(bins+(dr*0.5), star_sfr_lower[m], star_sfr_higher[m], alpha=0.3)
plt.legend()
plt.xlabel('R half *')
plt.ylabel('SFR surface density (Msun/yr /kpc^2)')
plt.xlim(0, 2)
plt.savefig(results_dir+'star_sfr_medians.png')
plt.clf()
for m in range(len(winds)):
	plt.errorbar(bins+(dr*0.5), star_sfr_mean[m], yerr=star_sfr_sigma[m], marker='.', markersize=4, linestyle='--', 
				label=wind_labels[m] +', '+str(no_gals[m])+' galaxies')
plt.legend()
plt.xlabel('R half *')
plt.ylabel('SFR surface density (Msun/yr /kpc^2)')
plt.xlim(0, 2)
plt.savefig(results_dir+'star_sfr_means.png')
plt.clf()


for m in range(len(winds)):
	plt.plot(bins+(dr*0.5), np.log10(star_ssfr_median[m]), marker='.', markersize=4, linestyle='--', label=wind_labels[m] +', '+str(no_gals[m])+' galaxies')
	plt.fill_between(bins+(dr*0.5), np.log10(star_ssfr_lower[m]), np.log10(star_ssfr_higher[m]), alpha=0.3)
plt.legend()
plt.xlabel('R half *')
plt.ylabel('log sSFR (yr^-1)')
plt.xlim(0, 2)
plt.ylim(-13, )
plt.savefig(results_dir+'star_ssfr_medians.png')
plt.clf()
for m in range(len(winds)):
	plt.errorbar(bins+(dr*0.5), np.log10(star_ssfr_mean[m]), yerr=star_ssfr_sigma[m], marker='.', markersize=4, linestyle='--', 
				label=wind_labels[m] +', '+str(no_gals[m])+' galaxies')
plt.legend()
plt.xlabel('R half *')
plt.ylabel('log sSFR (yr^-1)')
plt.xlim(0, 2)
plt.ylim(-13, )
plt.savefig(results_dir+'star_ssfr_means.png')
plt.clf()


for m in range(len(winds)):
	plt.plot(bins+(dr*0.5), gas_sfr_median[m], marker='.', markersize=4, linestyle='--', label=wind_labels[m] +', '+str(no_gals[m])+' galaxies')
	plt.fill_between(bins+(dr*0.5), gas_sfr_lower[m], gas_sfr_higher[m], alpha=0.3)
plt.legend()
plt.xlabel('R half *')
plt.ylabel('SFR surface density (Msun/yr /kpc^2)')
plt.xlim(0, 2)
plt.ylim(0, )
plt.savefig(results_dir+'gas_sfr_medians.png')
plt.clf()
for m in range(len(winds)):
	plt.errorbar(bins+(dr*0.5), gas_sfr_mean[m], yerr=gas_sfr_sigma[m], marker='.', markersize=4, linestyle='--', 
				label=wind_labels[m] +', '+str(no_gals[m])+' galaxies')
plt.legend()
plt.xlabel('R half *')
plt.ylabel('SFR surface density (Msun/yr /kpc^2)')
plt.xlim(0, 2)
plt.ylim(0, )
plt.savefig(results_dir+'gas_sfr_means.png')
plt.clf()


for m in range(len(winds)):
	plt.plot(bins+(dr*0.5), np.log10(gas_ssfr_median[m]), marker='.', markersize=4, linestyle='--', label=wind_labels[m] +', '+str(no_gals[m])+' galaxies')
	plt.fill_between(bins+(dr*0.5), np.log10(gas_ssfr_lower[m]), np.log10(gas_ssfr_higher[m]), alpha=0.3)
plt.legend()
plt.xlabel('R half *')
plt.ylabel('log sSFR (yr^-1)')
plt.xlim(0, 2)
plt.ylim(-13, )
plt.savefig(results_dir+'gas_ssfr_medians.png')
plt.clf()
for m in range(len(winds)):
	plt.errorbar(bins+(dr*0.5), np.log10(gas_ssfr_mean[m]), yerr=gas_ssfr_sigma[m], marker='.', markersize=4, linestyle='--', 
				label=wind_labels[m] +', '+str(no_gals[m])+' galaxies')
plt.legend()
plt.xlabel('R half *')
plt.ylabel('log sSFR (yr^-1)')
plt.xlim(0, 2)
plt.ylim(-13, )
plt.savefig(results_dir+'gas_ssfr_means.png')
plt.clf()


for m in range(len(winds)):
	plt.plot(bins+(dr*0.5), gas_h1_median[m], marker='.', markersize=4, linestyle='--', label=wind_labels[m] +', '+str(no_gals[m])+' galaxies')
	plt.fill_between(bins+(dr*0.5), gas_h1_lower[m], gas_h1_higher[m], alpha=0.3)
plt.legend()
plt.xlabel('R half *')
plt.ylabel('HI surface density (Msun /kpc^2)')
plt.xlim(0, 2)
plt.ylim(0, )
plt.savefig(results_dir+'h1_medians.png')
plt.clf()
for m in range(len(winds)):
	plt.errorbar(bins+(dr*0.5), gas_h1_mean[m], yerr=gas_h1_sigma[m], marker='.', markersize=4, linestyle='--', 
				label=wind_labels[m] +', '+str(no_gals[m])+' galaxies')
plt.legend()
plt.xlabel('R half *')
plt.ylabel('HI surface density (Msun /kpc^2)')
plt.xlim(0, 2)
plt.ylim(0, )
plt.savefig(results_dir+'h1_means.png')
plt.clf()


for m in range(len(winds)):
	plt.plot(bins+(dr*0.5), gas_h2_median[m], marker='.', markersize=4, linestyle='--', label=wind_labels[m] +', '+str(no_gals[m])+' galaxies')
	plt.fill_between(bins+(dr*0.5), gas_h2_lower[m], gas_h2_higher[m], alpha=0.3)
plt.legend()
plt.xlabel('R half *')
plt.ylabel('HII surface density (Msun /kpc^2)')
plt.xlim(0, 2)
plt.ylim(0, )
plt.savefig(results_dir+'h2_medians.png')
plt.clf()
for m in range(len(winds)):
	plt.errorbar(bins+(dr*0.5), gas_h2_mean[m], yerr=gas_h2_sigma[m], marker='.', markersize=4, linestyle='--', 
				label=wind_labels[m] +', '+str(no_gals[m])+' galaxies')
plt.legend()
plt.xlabel('R half *')
plt.ylabel('HII surface density (Msun /kpc^2)')
plt.xlim(0, 2)
plt.ylim(0, )
plt.savefig(results_dir+'h2_means.png')
plt.clf()