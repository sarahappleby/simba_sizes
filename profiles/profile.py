import h5py
import sys
import numpy as np 
import os
import gc
import matplotlib.pyplot as plt

sys.path.append('/home/sapple/tools/')
from projection import *

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
selection = 'green_valley'

n = 10 # number of bins
factor = 2.
dr = factor / n
rplot = np.arange(0., dr*n, dr) + (dr*0.5)

# for making profiles
mass_bins = [10., 10.5, 11.]
bin_labels = ['10.0 - 10.5', '10.5 - 11.0', '> 11.0']

model = sys.argv[1]
wind = sys.argv[2]
snap = sys.argv[3]
sample_file = sys.argv[4]
results_dir = sys.argv[5]

results_dir += '/'+model + '_' + snap + '/' + selection 
if rotate_galaxies:
	results_dir += '/rotated_faceon'
else:
	results_dir += '/random_orientation'
if bh_center:
	results_dir += '/bh_centered/'
else:
	results_dir += '/scom_centered/'
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

gal_sm = np.log10(gal_sm[gal_ids*gal_cent])
gal_ssfr = np.log10(gal_ssfr[gal_ids*gal_cent])
gal_rad = gal_rad[gal_ids*gal_cent]
gal_ids = np.arange(len(gal_ids))[gal_ids*gal_cent]

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

no_gals = np.zeros(3)

sm_median = np.zeros((3, n)); sm_lower = np.zeros((3, n)); sm_higher = np.zeros((3, n)); sm_mean = np.zeros((3, n)); sm_sigma = np.zeros((3, n))
star_sfr_median = np.zeros((3, n)); star_sfr_lower = np.zeros((3, n)); star_sfr_higher = np.zeros((3, n)); star_sfr_mean = np.zeros((3, n)); star_sfr_sigma = np.zeros((3, n))
star_ssfr_median = np.zeros((3, n)); star_ssfr_lower = np.zeros((3, n)); star_ssfr_higher = np.zeros((3, n)); star_ssfr_mean = np.zeros((3, n)); star_ssfr_sigma = np.zeros((3, n))

gas_sfr_median = np.zeros((3, n)); gas_sfr_lower = np.zeros((3, n)); gas_sfr_higher = np.zeros((3, n)); gas_sfr_mean = np.zeros((3, n)); gas_sfr_sigma = np.zeros((3, n))
gas_ssfr_median = np.zeros((3, n)); gas_ssfr_lower = np.zeros((3, n)); gas_ssfr_higher = np.zeros((3, n)); gas_ssfr_mean = np.zeros((3, n)); gas_ssfr_sigma = np.zeros((3, n))
gas_h1_median = np.zeros((3, n)); gas_h1_lower = np.zeros((3, n)); gas_h1_higher = np.zeros((3, n)); gas_h1_mean = np.zeros((3, n)); gas_h1_sigma = np.zeros((3, n))
gas_h2_median = np.zeros((3, n)); gas_h2_lower = np.zeros((3, n)); gas_h2_higher = np.zeros((3, n)); gas_h2_mean = np.zeros((3, n)); gas_h2_sigma = np.zeros((3, n))


for m in range(len(mass_bins)):
	print '\n'
	print 'Looking at mass bin ' + bin_labels[m]
	if m != 2:
		sm_mask = (gal_sm > mass_bins[m]) & (gal_sm < mass_bins[m+1])
	else:
		sm_mask = gal_sm > mass_bins[m]

	gal_sm_use = gal_sm[sm_mask]
	gal_rad_use = gal_rad[sm_mask]
	gal_ssfr_use = gal_ssfr[sm_mask]
	gal_ids_use = gal_ids[sm_mask]

	no_gals[m] = len(gal_ids_use)
	print str(no_gals[m]) + ' galaxies in bin'
	print '\n'

	use_star_m = np.zeros((len(gal_ids_use), n))
	use_star_sfr = np.zeros((len(gal_ids_use), n))
	use_star_ssfr = np.zeros((len(gal_ids_use), n))
	use_gas_sfr = np.zeros((len(gal_ids_use), n))
	use_gas_ssfr = np.zeros((len(gal_ids_use), n))
	use_gas_h1 = np.zeros((len(gal_ids_use), n))
	use_gas_h2 = np.zeros((len(gal_ids_use), n))


	for i in range(len(gal_ids_use)):

		print '\n'
		print 'Galaxy ' +str(gal_ids_use[i])
		slist = sim.galaxies[gal_ids_use[i]].slist
		glist = sim.galaxies[gal_ids_use[i]].glist
		rhalf = gal_rad_use[i]
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

		if rotate_galaxies:
			axis, angle = compute_rotation_to_vec(pos[:, 0], pos[:, 1], pos[:, 2], vel[:, 0], vel[:, 1], vel[:, 2], mass, vec)
			pos[:, 0], pos[:, 1], pos[:, 2] = rotate(pos[:, 0], pos[:, 1], pos[:, 2], axis, angle)
			r = np.linalg.norm(pos[:, [0, 1]], axis=1) / rhalf
		else:
			r = np.sqrt(pos[:, 0]**2 +pos[:, 1]**2) / rhalf

		plot_name = results_dir + 'images/gal_'+str(gal_ids_use[i]) + '_stars.png'
		make_image(pos[:, 0], pos[:, 1], np.log10(mass), plot_name, rhalf)

		use_star_m[i] = make_profile(n, dr, r, mass, rhalf)
		use_star_sfr[i] = make_profile(n, dr, r[ages_mask], mass[ages_mask], rhalf) * mass_loss / time

		plot_profile(rplot, np.log10(use_star_m[i]), results_dir+'profiles/sm_profile_gal_'+str(gal_ids_use[i])+'.png', 'M* surface density', title=title)
		plot_profile(rplot, use_star_sfr[i], results_dir+'profiles/star_sfr_profile_gal_'+str(gal_ids_use[i])+'.png', 'SFR surface density', title=title)

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
		Rotate to face on plane
		"""
		if rotate_galaxies:
			axis, angle = compute_rotation_to_vec(pos[:, 0], pos[:, 1], pos[:, 2], vel[:, 0], vel[:, 1], vel[:, 2], mass, vec)
			pos[:, 0], pos[:, 1], pos[:, 2] = rotate(pos[:, 0], pos[:, 1], pos[:, 2], axis, angle)
			r = np.linalg.norm(pos[:, [0, 1]], axis=1) / rhalf
		else:
			r  = np.sqrt(pos[:, 0]**2 +pos[:, 1]**2) / rhalf

		plot_name = results_dir + 'images/gal_'+str(gal_ids_use[i]) + '_gas.png'
		make_image(pos[:, 0], pos[:, 1], np.log10(mass), plot_name, rhalf)
		
		use_gas_sfr[i] = make_profile(n, dr, r, sfr, rhalf)
		use_gas_h1[i] = make_profile(n, dr, r, h1, rhalf)
		use_gas_h2[i] = make_profile(n, dr, r, h2, rhalf)

		use_gas_h1[i] /= np.sum(use_gas_h1[i])
		use_gas_h2[i] /= np.sum(use_gas_h2[i])

		plot_profile(rplot, use_gas_sfr[i], results_dir+'profiles/gas_sfr_profile_gal_'+str(gal_ids_use[i])+'.png', 'SFR surface density', title=title)
		plot_profile(rplot, use_gas_h1[i], results_dir+'profiles/gas_h1_profile_gal_'+str(gal_ids_use[i])+'.png', 'HI fraction surface density', title=title)
		plot_profile(rplot, use_gas_h2[i], results_dir+'profiles/gas_h2_profile_gal_'+str(gal_ids_use[i])+'.png', 'HII fraction surface density', title=title)

	with h5py.File(results_dir+'mask_'+str(m)+'_all_profiles.h5', 'a') as f:
		f.create_dataset('star_sfr', data=np.array(use_star_sfr))
		f.create_dataset('gas_sfr', data=np.array(use_gas_sfr))
		f.create_dataset('h1', data=np.array(use_gas_h1))
		f.create_dataset('h2', data=np.array(use_gas_h2))
		f.create_dataset('sm', data=np.array(use_star_m))
		f.create_dataset('gal_ids', data=np.array(gal_ids_use))
	
	del use_star_m, use_star_sfr, use_gas_sfr, use_gas_h1, use_gas_h2
	gc.collect()
	print '\n'
