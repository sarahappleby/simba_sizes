import h5py
import sys
import numpy as np 
import os

from readgadget import readsnap
import caesar
from astropy.cosmology import FlatLambdaCDM

from profile_methods import *

N = 5
age_min = 50.
age_max = 150.
mass_loss = 1.18
time = 100.e6
bh_center = True
rotate_galaxies = False

n = 10
factor = 2.
dr = factor / n
bins = np.arange(0., factor, dr)
#bins = np.zeros(n)
#bins[1:] = np.arange(dr*2., factor, dr)
area = np.array([np.pi*(i+dr)**2 - np.pi* i**2 for i in bins])

model = sys.argv[1]
wind = sys.argv[2]
snap = sys.argv[3]
sample_file = sys.argv[4]
results_dir = sys.argv[5]

if not os.path.exists(results_dir):
	os.makedirs(results_dir)
	os.makedirs(results_dir+'/images')
	os.makedirs(results_dir+'/profiles')

with h5py.File(sample_file, 'r') as f:
	try:
		gal_ids = f[model+'_'+snap].value
	except KeyError:
		print 'Need to identify galaxies; run gv_sample.py first'

data_dir = '/home/rad/data/'+model+'/'+wind+'/'
snapfile = data_dir+'snap_'+model+'_'+snap+'.hdf5'

sim =  caesar.load(data_dir+'Groups/'+model+'_'+snap+'.hdf5', LoadHalo=False)

h = sim.simulation.hubble_constant
redshift = sim.simulation.redshift
cosmo = FlatLambdaCDM(H0=100*h, Om0=sim.simulation.omega_matter, Ob0=sim.simulation.omega_baryon, Tcmb0=2.73)
thubble = cosmo.age(redshift).value # in Gyr

# load in the caesar galaxy data to make an initial cut of star forming galaxies
gal_cent = np.array([i.central for i in sim.galaxies])
gal_sm = np.array([i.masses['stellar'].in_units('Msun') for i in sim.galaxies])
gal_sfr = np.array([i.sfr.in_units('Msun/yr') for i in sim.galaxies])
gal_sfr[np.where(gal_sfr == 1.)[0]] = 0.
gal_ssfr = gal_sfr / gal_sm
gal_rad = np.array([i.radii['stellar_half_mass'].in_units('kpc') for i in sim.galaxies])


gal_sm_use = np.log10(gal_sm[gal_ids*gal_cent])
gal_ssfr_use = np.log10(gal_ssfr[gal_ids*gal_cent])
gal_rad_use = gal_rad[gal_ids*gal_cent]
gal_ids_use = np.arange(len(gal_ids))[gal_ids*gal_cent]

gal_nos = np.arange(N)


# load in star particle data with readsnap
star_pos = readsnap(snapfile, 'pos', 'star', suppress=1, units=1) / (h*(1.+redshift)) # in kpc
star_mass = readsnap(snapfile, 'mass', 'star', suppress=1, units=1) / h
star_vels = readsnap(snapfile, 'vel', 'star', suppress=1, units=0)
star_tform = readsnap(snapfile, 'age', 'star', suppress=1, units=1) # expansion times at time of formation

gas_pos = readsnap(snapfile, 'pos', 'gas', suppress=1, units=1) / (h*(1.+redshift)) # in kpc
gas_mass = readsnap(snapfile, 'mass', 'gas', suppress=1, units=1) / h
gas_vels = readsnap(snapfile, 'vel', 'gas', suppress=1, units=0)
gas_sfr = readsnap(snapfile, 'sfr', 'gas', suppress=1, units=1)
gas_h1 = readsnap(snapfile, 'NeutralHydrogenAbundance', 'gas', suppress=1, units=1)
gas_h2 = readsnap(snapfile, 'fh2', 'gas', suppress=1, units=1)

bh_pos = readsnap(snapfile, 'pos', 'bndry', suppress=1, units=1) / (h*(1.+redshift)) # in kpc

gal_sm_use = gal_sm_use[gal_nos]
gal_rad_use = gal_rad_use[gal_nos]
gal_ssfr_use = gal_ssfr_use[gal_nos]
gal_ids_use = gal_ids_use[gal_nos]

use_star_m = np.zeros((N, n))
use_star_sfr = np.zeros((N, n))
use_star_ssfr = np.zeros((N, n))
use_gas_sfr = np.zeros((N, n))
use_gas_ssfr = np.zeros((N, n))
use_gas_h1 = np.zeros((N, n))
use_gas_h2 = np.zeros((N, n))



for i in range(len(gal_ids_use)):
	print '\n'
	print 'Galaxy ' +str(gal_ids_use[i])
	slist = sim.galaxies[gal_ids_use[i]].slist
	glist = sim.galaxies[gal_ids_use[i]].glist
	r_max = factor*gal_rad_use[i]
	title = 'log M* = ' + str(round(gal_sm_use[i], 2)) + '; log(sSFR) = ' + format(round(gal_ssfr_use[i], 2))

	print str(len(glist)) + ' gas particles'
	print 'log sSFR: ' + format(round(gal_ssfr_use[i], 2))

	"""
	Find stellar ages:
	"""
	print 'Find stellar ages'
	# find ages of galaxy stars
	tform = star_tform[slist]
	ages = tage(cosmo, thubble, tform) * 1000. # get star ages in Myr
	ages_mask = (ages > age_min) & (ages < age_max)

	print 'Number of stars: ' +str(len(ages[ages_mask]))
	print 'Youngest star age: ' +str(np.min(ages))
	print 'Oldest star age: ' +str(np.max(ages))

	"""
	Get star particle data and correct for bh center or star com
	"""
	pos = star_pos[slist]
	vel = star_vels[slist]
	mass = star_mass[slist]
	if not bh_center:
		pos -= center_of_quantity(pos, mass)
	else:
		if len(sim.galaxies[gal_ids_use[i]].bhlist) > 0.:
			pos -= bh_pos[sim.galaxies[gal_ids_use[i]].bhlist[0]]
		else:
			pos -= center_of_quantity(pos, mass)
			print 'No black holes to center on, centering on stars'	
	vel -= center_of_quantity(vel, mass)
	"""
	Filter for particles within 2 half stellar mass radii
	"""
	r = np.linalg.norm(pos, axis=1)
	mass = mass[r < r_max ]
	vel = vel[r < r_max ]
	pos = pos[r < r_max ]
	ages_mask = np.array(ages_mask)[r < r_max]
	r = r[r < r_max]
	if len(r) == 0:
		print 'No particles within max radius - problem with center of mass'
		rad_mask[i] = False
		continue
	"""
	Rotate to face on plane
	"""
	if rotate_galaxies:
		axis, angle = compute_rotation_to_vec(pos[:, 0], pos[:, 1], pos[:, 2], vel[:, 0], vel[:, 1], vel[:, 2], mass, vec)
		pos[:, 0], pos[:, 1], pos[:, 2] = rotate(pos[:, 0], pos[:, 1], pos[:, 2], axis, angle)
		r = (np.linalg.norm(pos[:, [0, 1]], axis=1) / r_max) * 2.
	else:
		#r  = (r/ r_max)*2.
		r = r
	plot_name = results_dir + 'images/gal_'+str(gal_ids_use[i]) + '_stars.png'
	make_image(pos[:, 0], pos[:, 1], np.log10(mass), plot_name, r_max)

	use_star_m[i] = make_profile(n, dr, r, mass, 1)
	use_star_sfr[i] = make_profile(n, dr, r[ages_mask], mass[ages_mask], r_max) * mass_loss / time

	use_star_ssfr[i] = use_star_sfr[i] / use_star_m[i]
	use_star_ssfr[i][np.where(use_star_ssfr[i] == 0.)] = 1.e-20

	plot_profile(bins+(dr*0.5), use_star_m[i], results_dir+'profiles/sm_profile_gal_'+str(gal_ids_use[i])+'.png', 'M* surface density', title=title)
	plot_profile(bins+(dr*0.5), use_star_sfr[i], results_dir+'profiles/star_sfr_profile_gal_'+str(gal_ids_use[i])+'.png', 'SFR surface density', title=title)
	plot_profile(bins+(dr*0.5), np.log10(use_star_ssfr[i]), results_dir+'profiles/star_ssfr_profile_gal_'+str(gal_ids_use[i])+'.png', 'log sSFR', title=title, ylim=-13)

	"""
	For the gas particles:
	"""
	pos = gas_pos[glist]
	vel = gas_vels[glist]
	mass = gas_mass[glist]
	sfr = gas_sfr[glist]
	h1 = gas_h1[glist]
	h2 = gas_h2[glist]
	if not bh_center:
		pos -= center_of_quantity(star_pos[slist], star_mass[slist])
	else:
		if len(sim.galaxies[gal_ids_use[i]].bhlist) > 0.:
			pos -= bh_pos[sim.galaxies[gal_ids_use[i]].bhlist[0]]
		else:
			pos -= center_of_quantity(pos, mass)
			print 'No black holes to center on, centering on stars'	
	vel -= center_of_quantity(vel, mass)
	"""
	Filter for particles within 2 half stellar mass radii
	"""
	r = np.linalg.norm(pos, axis=1)
	mass = mass[r < r_max ]
	vel = vel[r < r_max ]
	pos = pos[r < r_max ]
	sfr = sfr[r < r_max]
	h1 = h1[r < r_max] * mass
	h2 = h2[r < r_max] * mass
	r = r[r < r_max]
	if len(r) != 0.:
		"""
		Rotate to face on plane
		"""
		if rotate_galaxies:
			axis, angle = compute_rotation_to_vec(pos[:, 0], pos[:, 1], pos[:, 2], vel[:, 0], vel[:, 1], vel[:, 2], mass, vec)
			pos[:, 0], pos[:, 1], pos[:, 2] = rotate(pos[:, 0], pos[:, 1], pos[:, 2], axis, angle)
			r = (np.linalg.norm(pos[:, [0, 1]], axis=1) / r_max) * 2.
		else:
			r  = (r/ r_max)*2.
		plot_name = results_dir + 'images/gal_'+str(gal_ids_use[i]) + '_gas.png'
		make_image(pos[:, 0], pos[:, 1], np.log10(mass), plot_name, r_max)
		
		use_gas_sfr[i] = make_profile(n, dr, r, sfr)
		use_gas_h1[i] = make_profile(n, dr, r, h1)
		use_gas_h2[i] = make_profile(n, dr, r, h2)

		use_gas_h1[i] /= np.sum(use_gas_h1[i])
		use_gas_h2[i] /= np.sum(use_gas_h2[i])

	use_gas_ssfr[i] = use_gas_sfr[i] / use_star_m[i]
	use_gas_ssfr[i][np.where(use_gas_ssfr[i] == 0.)] = 1.e-20

	plot_profile(bins+(dr*0.5), use_gas_sfr[i], results_dir+'profiles/gas_sfr_profile_gal_'+str(gal_ids_use[i])+'.png', 'SFR surface density', title=title)
	plot_profile(bins+(dr*0.5), np.log10(use_gas_ssfr[i]), results_dir+'profiles/gas_ssfr_profile_gal_'+str(gal_ids[i])+'.png', 'log sSFR', title=title, ylim=-13)
	plot_profile(bins+(dr*0.5), use_gas_h1[i], results_dir+'profiles/gas_h1_profile_gal_'+str(gal_ids_use[i])+'.png', 'HI fraction surface density', title=title)
	plot_profile(bins+(dr*0.5), use_gas_h2[i], results_dir+'profiles/gas_h2_profile_gal_'+str(gal_ids_use[i])+'.png', 'HII fraction surface density', title=title)


sm_median = np.percentile(use_star_m, '50', axis=0)
sm_lower = np.percentile(use_star_m, '25', axis=0)
sm_higher = np.percentile(use_star_m, '75', axis=0)
sm_mean = np.mean(use_star_m, axis=0)
sm_sigma = np.std(use_star_m, axis=0)

star_sfr_median = np.nanpercentile(use_star_sfr, '50', axis=0)
star_sfr_lower = np.nanpercentile(use_star_sfr, '25', axis=0)
star_sfr_higher = np.nanpercentile(use_star_sfr, '75', axis=0)
star_sfr_mean = np.mean(use_star_sfr, axis=0)
star_sfr_sigma = np.std(use_star_sfr, axis=0)

# in manga paper, ssfr profiles are ratio of median sfr and median sm

star_ssfr_median = star_sfr_median / sm_median
star_ssfr_lower = star_sfr_lower / sm_lower
star_ssfr_higher= star_sfr_higher / sm_higher
star_ssfr_mean = star_sfr_mean/ sm_mean
star_ssfr_sigma = star_ssfr_mean * ((star_sfr_sigma / star_sfr_mean)**2 + (sm_sigma / sm_mean)**2)**0.5

star_ssfr_median[star_ssfr_median == 0.] = 1.e-20
star_ssfr_lower[star_ssfr_lower == 0.] = 1.e-20
star_ssfr_higher[star_ssfr_higher == 0.] = 1.e-20
star_ssfr_mean[star_ssfr_mean == 0.] = 1.e-20
star_ssfr_sigma /= (star_ssfr_mean*np.log(10.))


gas_sfr_median = np.nanpercentile(use_gas_sfr, '50', axis=0)
gas_sfr_lower = np.nanpercentile(use_gas_sfr, '25', axis=0)
gas_sfr_higher = np.nanpercentile(use_gas_sfr, '75', axis=0)
gas_sfr_mean = np.mean(use_gas_sfr, axis=0)
gas_sfr_sigma = np.std(use_gas_sfr, axis=0)

gas_ssfr_median = gas_sfr_median / sm_median
gas_ssfr_lower = gas_sfr_lower / sm_lower
gas_ssfr_higher = gas_sfr_higher / sm_higher
gas_ssfr_mean = gas_sfr_mean / sm_mean
gas_ssfr_sigma = gas_ssfr_mean * ((gas_sfr_sigma / gas_sfr_mean)**2 + (sm_sigma / sm_mean)**2)**0.5

gas_ssfr_median[gas_ssfr_median == 0.] = 1.e-20
gas_ssfr_lower[gas_ssfr_lower == 0.] = 1.e-20
gas_ssfr_higher[gas_ssfr_higher == 0.] = 1.e-20
gas_ssfr_mean[gas_ssfr_mean == 0.] = 1.e-20
gas_ssfr_sigma /= (gas_ssfr_mean*np.log(10.))

gas_h1_median = np.nanpercentile(use_gas_h1, '50', axis=0)
gas_h1_lower = np.nanpercentile(use_gas_h1, '25', axis=0)
gas_h1_higher = np.nanpercentile(use_gas_h1, '75', axis=0)
gas_h1_mean = np.nanmean(use_gas_h1, axis=0)
gas_h1_sigma = np.nanstd(use_gas_h1, axis=0)

gas_h2_median = np.nanpercentile(use_gas_h2, '50', axis=0)
gas_h2_lower = np.nanpercentile(use_gas_h2, '25', axis=0)
gas_h2_higher = np.nanpercentile(use_gas_h2, '75', axis=0)
gas_h2_mean = np.nanmean(use_gas_h2, axis=0)
gas_h2_sigma = np.nanstd(use_gas_h2, axis=0)

plt.semilogy(bins+(dr*0.5), sm_median, marker='.', markersize=4, linestyle='--')
plt.fill_between(bins+(dr*0.5), sm_lower, sm_higher, alpha=0.3)
plt.xlabel('R half *')
plt.ylabel('M* surface density (Msun/R half bin^2)')
plt.xlim(0, 2)
plt.savefig(results_dir+'sm_medians.png')
plt.clf()
plt.errorbar(bins+(dr*0.5), sm_mean, yerr=sm_sigma, marker='.', markersize=4, linestyle='--')
plt.yscale('log')
plt.xlabel('R half *')
plt.ylabel('M* surface density (Msun/R half bin^2)')
plt.xlim(0, 2)
plt.savefig(results_dir+'sm_means.png')
plt.clf()


plt.plot(bins+(dr*0.5), star_sfr_median, marker='.', markersize=4, linestyle='--')
plt.fill_between(bins+(dr*0.5), star_sfr_lower, star_sfr_higher, alpha=0.3)
plt.xlabel('R half *')
plt.ylabel('SFR surface density (Msun/yr /R half bin^2)')
plt.xlim(0, 2)
plt.savefig(results_dir+'star_sfr_medians.png')
plt.clf()
plt.errorbar(bins+(dr*0.5), star_sfr_mean, yerr=star_sfr_sigma, marker='.', markersize=4, linestyle='--')
plt.xlabel('R half *')
plt.ylabel('SFR surface density (Msun/yr /R half bin^2)')
plt.xlim(0, 2)
plt.savefig(results_dir+'star_sfr_means.png')
plt.clf()


plt.plot(bins+(dr*0.5), np.log10(star_ssfr_median), marker='.', markersize=4, linestyle='--')
plt.fill_between(bins+(dr*0.5), np.log10(star_ssfr_lower), np.log10(star_ssfr_higher), alpha=0.3)
plt.xlabel('R half *')
plt.ylabel('log sSFR (yr^-1)')
plt.xlim(0, 2)
plt.ylim(-13, )
plt.savefig(results_dir+'star_ssfr_medians.png')
plt.clf()
plt.errorbar(bins+(dr*0.5), np.log10(star_ssfr_mean), yerr=star_ssfr_sigma, marker='.', markersize=4, linestyle='--')
plt.xlabel('R half *')
plt.ylabel('log sSFR (yr^-1)')
plt.xlim(0, 2)
plt.ylim(-13, )
plt.savefig(results_dir+'star_ssfr_means.png')
plt.clf()


plt.plot(bins+(dr*0.5), gas_sfr_median, marker='.', markersize=4, linestyle='--')
plt.fill_between(bins+(dr*0.5), gas_sfr_lower, gas_sfr_higher, alpha=0.3)
plt.xlabel('R half *')
plt.ylabel('SFR surface density (Msun/yr /R half bin^2)')
plt.xlim(0, 2)
plt.ylim(0, )
plt.savefig(results_dir+'gas_sfr_medians.png')
plt.clf()
plt.errorbar(bins+(dr*0.5), gas_sfr_mean, yerr=gas_sfr_sigma, marker='.', markersize=4, linestyle='--')
plt.xlabel('R half *')
plt.ylabel('SFR surface density (Msun/yr /R half bin^2)')
plt.xlim(0, 2)
plt.ylim(0, )
plt.savefig(results_dir+'gas_sfr_means.png')
plt.clf()


plt.plot(bins+(dr*0.5), np.log10(gas_ssfr_median), marker='.', markersize=4, linestyle='--')
plt.fill_between(bins+(dr*0.5), np.log10(gas_ssfr_lower), np.log10(gas_ssfr_higher), alpha=0.3)
plt.xlabel('R half *')
plt.ylabel('log sSFR (yr^-1)')
plt.xlim(0, 2)
plt.ylim(-13, )
plt.savefig(results_dir+'gas_ssfr_medians.png')
plt.clf()
plt.errorbar(bins+(dr*0.5), np.log10(gas_ssfr_mean), yerr=gas_ssfr_sigma, marker='.', markersize=4, linestyle='--')
plt.xlabel('R half *')
plt.ylabel('log sSFR (yr^-1)')
plt.xlim(0, 2)
plt.ylim(-13, )
plt.savefig(results_dir+'gas_ssfr_means.png')
plt.clf()


plt.plot(bins+(dr*0.5), gas_h1_median, marker='.', markersize=4, linestyle='--')
plt.fill_between(bins+(dr*0.5), gas_h1_lower, gas_h1_higher, alpha=0.3)
plt.xlabel('R half *')
plt.ylabel('HI surface density (Msun /R half bin^2)')
plt.xlim(0, 2)
plt.ylim(0, )
plt.savefig(results_dir+'h1_medians.png')
plt.clf()
plt.errorbar(bins+(dr*0.5), gas_h1_mean, yerr=gas_h1_sigma, marker='.', markersize=4, linestyle='--')
plt.xlabel('R half *')
plt.ylabel('HI surface density (Msun /R half bin^2)')
plt.xlim(0, 2)
plt.ylim(0, )
plt.savefig(results_dir+'h1_means.png')
plt.clf()


plt.plot(bins+(dr*0.5), gas_h2_median, marker='.', markersize=4, linestyle='--')
plt.fill_between(bins+(dr*0.5), gas_h2_lower, gas_h2_higher, alpha=0.3)
plt.xlabel('R half *')
plt.ylabel('HII surface density (Msun /R half bin^2)')
plt.xlim(0, 2)
plt.ylim(0, )
plt.savefig(results_dir+'h2_medians.png')
plt.clf()
plt.errorbar(bins+(dr*0.5), gas_h2_mean, yerr=gas_h2_sigma, marker='.', markersize=4, linestyle='--')
plt.xlabel('R half *')
plt.ylabel('HII surface density (Msun /R half bin^2)')
plt.xlim(0, 2)
plt.ylim(0, )
plt.savefig(results_dir+'h2_means.png')
plt.clf()