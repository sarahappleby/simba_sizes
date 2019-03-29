import sys
sys.path.append('/home/sapple/tools/')
from projection import *

import numpy as np 
from readgadget import readsnap
import caesar
import pygad as pg
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
import h5py
import os
import gc

from profile_methods import *

model = 'm50n512'
wind = 's50j7k'
snap = '151'

results_dir = '/home/sapple/simba_sizes/profiles/'+model+'/snap_'+snap+'/'
data_dir = '/home/rad/data/'+model+'/'+wind+'/'
snapfile = data_dir+'snap_'+model+'_'+snap+'.hdf5'

s = pg.Snap(snapfile)
sim =  caesar.load(data_dir+'Groups/'+model+'_'+snap+'.hdf5', LoadHalo=False)

h = sim.simulation.hubble_constant
redshift = sim.simulation.redshift
cosmo = FlatLambdaCDM(H0=100*h, Om0=sim.simulation.omega_matter, Ob0=sim.simulation.omega_baryon, Tcmb0=2.73)
thubble = cosmo.age(redshift).value # in Gyr

# for sample of galaxies and stars
ssfr_min = 3.e-10
sm_min = 5.e9
age_min = 50.
age_max = 150.
mass_loss = 1.18
time = 100.e6
do_ages = True
nstar_min = 100
bh_center = True

# for making profiles
DR = 1./(h*(1+redshift))
factor = 2.
vec = np.array([0, 0, 1]) # face-on projection to collapse the z direction
bins = np.arange(0., 2., 0.2)

# load in the caesar galaxy data to make an initial cut of star forming galaxies
gal_cent = np.array([i.central for i in sim.galaxies])
gal_sm = np.array([i.masses['stellar'].in_units('Msun') for i in sim.galaxies])
gal_bm = np.array([i.masses['bh'].in_units('Msun') for i in sim.galaxies])
gal_pos = np.array([i.pos.in_units('kpc') for i in sim.galaxies])
#gal_rad = np.array([i.radii['total_half_mass'].in_units('kpc') for i in sim.galaxies])
gal_sfr = np.array([i.sfr.in_units('Msun/yr') for i in sim.galaxies])
gal_ssfr = gal_sfr / gal_sm

sm_mask = gal_sm > sm_min
ssfr_mask = gal_ssfr > ssfr_min
gal_ids = np.arange(0, len(sim.galaxies))[ssfr_mask*sm_mask*gal_cent]

print 'Identified ' + str(len(gal_ids)) + ' galaxies'

rad_file = '/home/sapple/simba_sizes/sizes/data/'+model+'_'+wind+'_'+snap+'_data.h5'
if os.path.isfile(rad_file):
	with h5py.File(rad_file, 'r') as f:
		gal_rad = f['halfmass'].value
else:
	gal_rad = np.array([i.radii['stellar_half_mass'].in_units('kpc') for i in sim.galaxies])

# load in star particle data with readsnap
star_pos = readsnap(snapfile, 'pos', 'star', suppress=1, units=1) / (h*(1.+redshift)) # in kpc
star_mass = readsnap(snapfile, 'mass', 'star', suppress=1, units=1) / h
star_vels = readsnap(snapfile, 'vel', 'star', suppress=1, units=0)
star_tform = readsnap(snapfile, 'age', 'star', suppress=1, units=1) # expansion times at time of formation

gas_pos = readsnap(snapfile, 'pos', 'gas', suppress=1, units=1) / (h*(1.+redshift)) # in kpc
gas_mass = readsnap(snapfile, 'mass', 'gas', suppress=1, units=1) / h
gas_vels = readsnap(snapfile, 'vel', 'gas', suppress=1, units=0)
gas_sfr = readsnap(snapfile, 'sfr', 'gas', suppress=1, units=1)
abundance_h1 = s.gas['NeutralHydrogenAbundance']
abundance_h2 = s.gas['fh2']

bh_all_pos = s.bh['pos'].in_units_of('kpc')

gas_sfr_profiles = np.zeros((len(gal_ids), 10))
gas_ssfr_profiles = np.zeros((len(gal_ids), 10))
h1_profiles = np.zeros((len(gal_ids), 10))
h2_profiles = np.zeros((len(gal_ids), 10))
star_sfr_profiles = np.zeros((len(gal_ids), 10))
star_ssfr_profiles = np.zeros((len(gal_ids), 10))
sm_profiles = np.zeros((len(gal_ids), 10))
sm_frac_profiles = np.zeros((len(gal_ids), 10))

len_mask = np.array([True for i in range(len(gal_ids))])
rad_mask = np.array([True for i in range(len(gal_ids))])

for i in range(len(gal_ids)):

	print 'Galaxy ' +str(gal_ids[i])
	slist = sim.galaxies[gal_ids[i]].slist
	glist = sim.galaxies[gal_ids[i]].glist
	r_max = factor*gal_rad[gal_ids[i]]

	# for the profile
	NR = int(round(r_max / DR))
	if NR < 5:
		rad_mask[i] = False
		print 'Galaxy ' + str(gal_ids[i]) + ' too small (below 5 bins)'
		continue
	
	"""
	For the star particles:
	"""

	if do_ages:
		print 'Find stellar ages'
		# find ages of galaxy stars
		tform = star_tform[slist]
		ages = tage(cosmo, thubble, tform) * 1000. # get star ages in Myr
		ages_mask = (ages > age_min) & (ages < age_max)

		print 'Number of stars: ' +str(len(ages[ages_mask]))
		print 'Youngest star age: ' +str(np.min(ages))
		print 'Oldest star age: ' +str(np.max(ages))

		# filter for galaxies with more than 256 new star particles for the profiles
		if len(ages[ages_mask]) < nstar_min:
			len_mask[i] = False
	else:
		ages_mask = [True for j in range(len(slist))]

	pos = star_pos[slist]
	vel = star_vels[slist]
	mass = star_mass[slist]

	if not bh_center:
		pos -= center_of_quantity(pos, mass)
	else:
		pos -= bh_all_pos[sim.galaxies[gal_ids[i]].bhlist[0]]
	vel -= center_of_quantity(vel, mass)

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

	axis, angle = compute_rotation_to_vec(pos[:, 0], pos[:, 1], pos[:, 2], vel[:, 0], vel[:, 1], vel[:, 2], mass, vec)
	pos[:, 0], pos[:, 1], pos[:, 2] = rotate(pos[:, 0], pos[:, 1], pos[:, 2], axis, angle)
	plot_name = results_dir + 'images/gal_'+str(gal_ids[i]) + '_stars.png'
	make_image(pos[:, 0], pos[:, 1], np.log10(mass), plot_name, r_max)
	r = np.linalg.norm(pos[:, [0, 1]], axis=1)

	sm_profile = make_profile(NR, DR, r, mass)
	young_profile = make_profile(NR, DR, r[ages_mask], mass[ages_mask])
	# divide by 50Myr to get SFR in last 50Myr with 18% instantaneous mass loss
	star_sfr_profile = young_profile*mass_loss / time
	star_ssfr_profile = star_sfr_profile / sm_profile
	sm_frac = sm_profile / np.sum(mass)

	"""
	For the gas particles:
	"""
	pos = gas_pos[glist]
	vel = gas_vels[glist]
	mass = gas_mass[glist]
	sfr = gas_sfr[glist]
	h1 = abundance_h1[glist]
	h2 = abundance_h2[glist]

	if not bh_center:
		pos -= center_of_quantity(star_pos[slist], star_mass[slist])
	else:
		pos -= bh_all_pos[sim.galaxies[gal_ids[i]].bhlist[0]]
	vel -= center_of_quantity(vel, mass)

	r = np.linalg.norm(pos, axis=1)
	mass = mass[r < r_max ]
	vel = vel[r < r_max ]
	pos = pos[r < r_max ]
	sfr = sfr[r < r_max]
	h1 = h1[r < r_max]
	h2 = h2[r < r_max]
	r = r[r < r_max]

	axis, angle = compute_rotation_to_vec(pos[:, 0], pos[:, 1], pos[:, 2], vel[:, 0], vel[:, 1], vel[:, 2], mass, vec)
	pos[:, 0], pos[:, 1], pos[:, 2] = rotate(pos[:, 0], pos[:, 1], pos[:, 2], axis, angle)
	plot_name = results_dir + 'images/gal_'+str(gal_ids[i]) + '_gas.png'
	make_image(pos[:, 0], pos[:, 1], np.log10(mass), plot_name, r_max)
	r = np.linalg.norm(pos[:, [0, 1]], axis=1)

	gas_sfr_profile = make_profile(NR, DR, r, sfr)
	gas_ssfr_profile = gas_sfr_profile / sm_profile
	h1_profile = make_profile(NR, DR, r, h1*mass)
	h2_profile = make_profile(NR, DR, r, h2*mass)


	"""
	Plotting and binning:
	"""

	r_plot = np.arange(0.5*DR, (NR+0.5)*DR, DR) # bin centers
	if len(r_plot) > NR:
		r_plot = r_plot[1:]

	plt.semilogy(r_plot,sm_profile, linestyle='--', marker='.')
	plt.xlabel('R (kpc)')
	plt.ylabel('M* surface density (Msun/kpc^2)')
	plt.xlim(0, )
	plt.savefig(results_dir+'profiles/sm_profile_gal_'+str(gal_ids[i])+'.png')
	plt.clf()

	plt.semilogy(r_plot,sm_frac, linestyle='--', marker='.')
	plt.xlabel('R (kpc)')
	plt.ylabel('M* fraction')
	plt.xlim(0, )
	plt.savefig(results_dir+'profiles/sm_frac_profile_gal_'+str(gal_ids[i])+'.png')
	plt.clf()	

	plt.plot(r_plot,star_sfr_profile, linestyle='--', marker='.')
	plt.xlabel('R (kpc)')
	plt.ylabel('SFR surface density (Msun/yr kpc^2)')
	plt.xlim(0, )
	plt.savefig(results_dir+'profiles/star_sfr_profile_gal_'+str(gal_ids[i])+'.png')
	plt.clf()

	plt.plot(r_plot,np.log10(star_ssfr_profile), linestyle='--', marker='.')
	plt.xlabel('R (kpc)')
	plt.ylabel('log sSFR (yr^-1')
	plt.xlim(0, )
	plt.savefig(results_dir+'profiles/star_ssfr_profile_gal_'+str(gal_ids[i])+'.png')
	plt.clf()	

	plt.plot(r_plot,gas_sfr_profile, linestyle='--', marker='.')
	plt.xlabel('R (kpc)')
	plt.ylabel('SFR surface density (Msun/yr kpc^2)')
	plt.xlim(0, )
	plt.savefig(results_dir+'profiles/gas_sfr_profile_gal_'+str(gal_ids[i])+'.png')
	plt.clf()

	plt.plot(r_plot,np.log10(gas_ssfr_profile), linestyle='--', marker='.')
	plt.xlabel('R (kpc)')
	plt.ylabel('log sSFR (yr^-1)')
	plt.xlim(0, )
	plt.savefig(results_dir+'profiles/gas_ssfr_profile_gal_'+str(gal_ids[i])+'.png')
	plt.clf()

	plt.plot(r_plot,h1_profile, linestyle='--', marker='.')
	plt.xlabel('R (kpc)')
	plt.ylabel('HI surface density (Msun / kpc^2)')
	plt.xlim(0, )
	plt.savefig(results_dir+'profiles/h1_profile_gal_'+str(gal_ids[i])+'.png')
	plt.clf()

	plt.plot(r_plot,h2_profile, linestyle='--', marker='.')
	plt.xlabel('R (kpc)')
	plt.ylabel('HII surface density (Msun / kpc^2)')
	plt.xlim(0, )
	plt.savefig(results_dir+'profiles/h2_profile_gal_'+str(gal_ids[i])+'.png')
	plt.clf()

	x = np.arange(0., r_max, r_max/NR) / r_max*factor
	digitized = np.digitize(x, bins)
	sm_bin = [sm_profile[digitized == j] for j in range(1, len(bins) +1)]
	star_sfr_bin = [star_sfr_profile[digitized == j] for j in range(1, len(bins) +1)]
	star_ssfr_bin = [star_ssfr_profile[digitized == j] for j in range(1, len(bins) +1)]
	sm_frac_bin = [sm_frac[digitized == j] for j in range(1, len(bins) +1)]
	gas_sfr_bin = [gas_sfr_profile[digitized == j] for j in range(1, len(bins) + 1)]
	gas_ssfr_bin = [gas_ssfr_profile[digitized == j] for j in range(1, len(bins) + 1)]
	h1_bin = [h1_profile[digitized == j] for j in range(1, len(bins) + 1)]
	h2_bin = [h2_profile[digitized == j] for j in range(1, len(bins) + 1)]

	star_sfr_profiles[i] = np.array([np.sum(j) for j in star_sfr_bin])
	star_ssfr_profiles[i] = np.array([np.sum(j) for j in star_ssfr_bin])
	sm_profiles[i] = np.array([np.sum(j) for j in sm_bin])
	sm_frac_profiles[i] = np.array([np.sum(j) for j in sm_frac_bin])
	gas_sfr_profiles[i] = np.array([np.sum(j) for j in gas_sfr_bin])
	gas_ssfr_profiles[i] = np.array([np.sum(j) for j in gas_ssfr_bin])
	h1_profiles[i] = np.array([np.sum(j) for j in h1_bin])	
	h2_profiles[i] = np.array([np.sum(j) for j in h2_bin])

"""

def filter_zeros(array):
	array = np.transpose(array)
	array = [i[np.nonzero(i)] for i in array]
	return array

sm_profiles = filter_zeros(sm_profiles[rad_mask])
sm_frac_profiles = filter_zeros(sm_frac_profiles[rad_mask])
gas_sfr_profiles = filter_zeros(gas_sfr_profiles[rad_mask])
gas_ssfr_profiles = filter_zeros(gas_ssfr_profiles[rad_mask])
h1_profiles = filter_zeros(h1_profiles[rad_mask])
h2_profiles = filter_zeros(h2_profiles[rad_mask])

gal_ids_sm = gal_ids[rad_mask]
gal_ids_sfr = gal_ids[len_mask*rad_mask]

star_sfr_profiles = star_sfr_profiles[len_mask*rad_mask]
star_ssfr_profiles = star_ssfr_profiles[len_mask*rad_mask]

if len(star_sfr_profiles) != 0:
	star_sfr_lower = [np.percentile(i, 25) for i in star_sfr_profiles]
	star_sfr_higher = [np.percentile(i, 75) for i in star_sfr_profiles]
	star_sfr_median = [np.percentile(i, 50) for i in star_sfr_profiles]

	star_ssfr_lower = [np.percentile(i, 25) for i in np.log10(star_ssfr_profiles)]
	star_ssfr_higher = [np.percentile(i, 75) for i in np.log10(star_ssfr_profiles)]
	star_ssfr_median = [np.percentile(i, 50) for i in np.log10(star_ssfr_profiles)]

	plt.plot(bins+0.1, star_sfr_median, marker='.', markersize=2)
	plt.xlabel('R half *')
	plt.ylabel('SFR surface density (Msun /yr kpc^2)')
	plt.xlim(0, )
	plt.fill_between(bins+0.1, star_sfr_lower, star_sfr_higher, facecolor='blue', alpha=0.3)
	plt.title(str(len(gal_ids_sfr))+' galaxies')
	plt.savefig(results_dir+'star_sfr_profile.png')
	plt.clf()

	plt.plot(bins+0.1, star_ssfr_median, marker='.', markersize=2)
	plt.xlabel('R half *')
	plt.ylabel('log sSFR (yr^-1)')
	plt.xlim(0, )
	plt.fill_between(bins+0.1, star_ssfr_lower, star_ssfr_higher, facecolor='blue', alpha=0.3)
	plt.title(str(len(gal_ids_sfr))+' galaxies')
	plt.savefig(results_dir+'star_ssfr_profile.png')
	plt.clf()

sm_lower = [np.percentile(i, 25) for i in sm_profiles]
sm_higher = [np.percentile(i, 75) for i in sm_profiles]
sm_median = [np.percentile(i, 50) for i in sm_profiles]

sm_frac_lower = [np.percentile(i, 25) for i in sm_frac_profiles]
sm_frac_higher = [np.percentile(i, 75) for i in sm_frac_profiles]
sm_frac_median = [np.percentile(i, 50) for i in sm_frac_profiles]

gas_sfr_lower = [np.percentile(i, 25) for i in gas_sfr_profiles]
gas_sfr_higher = [np.percentile(i, 75) for i in gas_sfr_profiles]
gas_sfr_median = [np.percentile(i, 50) for i in gas_sfr_profiles]

gas_ssfr_lower = [np.percentile(i, 25) for i in np.log10(gas_ssfr_profiles)]
gas_ssfr_higher = [np.percentile(i, 75) for i in np.log10(gas_ssfr_profiles)]
gas_ssfr_median = [np.percentile(i, 50) for i in np.log10(gas_ssfr_profiles)]

h1_lower = [np.percentile(i, 25) for i in h1_profiles]
h1_higher = [np.percentile(i, 75) for i in h1_profiles]
h1_median = [np.percentile(i, 50) for i in h1_profiles]

h2_lower = [np.percentile(i, 25) for i in h2_profiles]
h2_higher = [np.percentile(i, 75) for i in h2_profiles]
h2_median = [np.percentile(i, 50) for i in h2_profiles]

plt.semilogy(bins+0.1, sm_median, marker='.', markersize=2)
plt.xlabel('R half *')
plt.ylabel('M* surface density (Msun/kpc^2)')
plt.xlim(0, )
plt.fill_between(bins+0.1, sm_lower, sm_higher, facecolor='blue', alpha=0.3)
plt.title(str(len(gal_ids_sm))+' galaxies')
plt.savefig(results_dir+'sm_profile.png')
plt.clf()

plt.semilogy(bins+0.1, sm_frac_median, marker='.', markersize=2)
plt.xlabel('R half *')
plt.ylabel('M* surface density fraction (kpc^-2)')
plt.xlim(0, )
plt.fill_between(bins+0.1, sm_frac_lower, sm_frac_higher, facecolor='blue', alpha=0.3)
plt.title(str(len(gal_ids_sm))+' galaxies')
plt.savefig(results_dir+'sm_frac_profile.png')
plt.clf()

plt.plot(bins+0.1, gas_sfr_median, marker='.', markersize=2)
plt.xlabel('R half *')
plt.ylabel('SFR surface density (Msun /yr kpc^2)')
plt.xlim(0, )
plt.fill_between(bins+0.1, gas_sfr_lower, gas_sfr_higher, facecolor='blue', alpha=0.3)
plt.title(str(len(gal_ids_sm))+' galaxies')
plt.savefig(results_dir+'gas_sfr_profile.png')
plt.clf()

plt.plot(bins+0.1, gas_ssfr_median, marker='.', markersize=2)
plt.xlabel('R half *')
plt.ylabel('log sSFR (yr^-1)')
plt.xlim(0, )
plt.fill_between(bins+0.1, gas_ssfr_lower, gas_ssfr_higher, facecolor='blue', alpha=0.3)
plt.title(str(len(gal_ids_sm))+' galaxies')
plt.savefig(results_dir+'gas_ssfr_profile.png')
plt.clf()

plt.plot(bins+0.1, h1_median, marker='.', markersize=2)
plt.xlabel('R half *')
plt.ylabel('HI surface density (Msun kpc^2)')
plt.xlim(0, )
plt.fill_between(bins+0.1, h1_lower, h1_higher, facecolor='blue', alpha=0.3)
plt.title(str(len(gal_ids_sm))+' galaxies')
plt.savefig(results_dir+'h1_profile.png')
plt.clf()

plt.plot(bins+0.1, h2_median, marker='.', markersize=2)
plt.xlabel('R half *')
plt.ylabel('HII surface density (Msun kpc^2)')
plt.xlim(0, )
plt.fill_between(bins+0.1, h2_lower, h2_higher, facecolor='blue', alpha=0.3)
plt.title(str(len(gal_ids_sm))+' galaxies')
plt.savefig(results_dir+'h2_profile.png')
plt.clf()

with h5py.File(results_dir+model+'_'+snap+'_profiles.h5', 'a') as f:
	f.create_dataset('star_sfr_profile', data=np.array(star_sfr_profiles))
	f.create_dataset('star_ssfr_profile', data=np.array(star_ssfr_profiles))
	f.create_dataset('gas_sfr_profile', data=np.array(gas_sfr_profiles))
	f.create_dataset('gas_ssfr_profile', data=np.array(gas_ssfr_profiles))
	f.create_dataset('h1_profile', data=np.array(h1_profiles))
	f.create_dataset('h2_profile', data=np.array(h2_profiles))
	f.create_dataset('sm_profile', data=np.array(sm_profiles))
	f.create_dataset('gal_ids', data=np.array(gal_ids_sm))
	f.create_dataset('gal_ids_young_stars', data=np.array(gal_ids_sfr))
	f.create_dataset('sm_frac_profiles', data=np.array(sm_frac_profiles))
"""