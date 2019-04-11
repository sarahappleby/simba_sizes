import sys
sys.path.append('/home/sapple/tools/')
from projection import *

import numpy as np 
from readgadget import readsnap
import caesar
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
import h5py
import os
import gc

from profile_methods import *

model = sys.argv[1]
wind = sys.argv[2]
snap = sys.argv[3]

results_dir = '/home/sapple/simba_sizes/profiles/suppressed_galaxies/'+model+'_'+snap + '/'
if not os.path.exists(results_dir):
	os.makedirs(results_dir)
	os.makedirs(results_dir+'images/')
	os.makedirs(results_dir+'profiles/')
data_dir = '/home/rad/data/'+model+'/'+wind+'/'
snapfile = data_dir+'snap_'+model+'_'+snap+'.hdf5'

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
gas_h1 = readsnap(snapfile, 'NeutralHydrogenAbundance', 'gas', suppress=1, units=1)
gas_h2 = readsnap(snapfile, 'fh2', 'gas', suppress=1, units=1)

bh_pos = readsnap(snapfile, 'pos', 'bndry', suppress=1, units=1) / (h*(1.+redshift)) # in kpc

suppress_mask = np.array([False for i in range(len(gal_ids))])
rad_mask = np.array([True for i in range(len(gal_ids))])
len_mask = np.array([True for i in range(len(gal_ids))])

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
		pos -= bh_pos[sim.galaxies[gal_ids[i]].bhlist[0]]
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
	keep_pos = pos.copy(); keep_mass = mass.copy()
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
	h1 = gas_h1[glist]
	h2 = gas_h2[glist]

	if not bh_center:
		pos -= center_of_quantity(star_pos[slist], star_mass[slist])
	else:
		pos -= bh_pos[sim.galaxies[gal_ids[i]].bhlist[0]]
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
	r = np.linalg.norm(pos[:, [0, 1]], axis=1)

	gas_sfr_profile = make_profile(NR, DR, r, sfr)
	gas_ssfr_profile = gas_sfr_profile / sm_profile

	h1_profile = make_profile(NR, DR, r, h1*mass)
	h2_profile = make_profile(NR, DR, r, h2*mass)

	"""
	Identify galaxies with central suppresion 
	"""
	if (gas_ssfr_profile[0] < np.max(gas_ssfr_profile)) & (gas_ssfr_profile[1] < np.max(gas_ssfr_profile)):
		print "Galaxy "+ str(gal_ids[i]) + ' is centrally suppressed'
		suppress_mask[i] = True
	else:
		continue

	plot_name = results_dir + 'images/gal_'+str(gal_ids[i]) + '_gas.png'
	make_image(pos[:, 0], pos[:, 1], np.log10(mass), plot_name, r_max)

	plot_name = results_dir + 'images/gal_'+str(gal_ids[i]) + '_stars.png'
	make_image(keep_pos[:, 0], keep_pos[:, 1], np.log10(keep_mass), plot_name, r_max)

	"""
	Plotting and binning:
	"""

	r_plot = np.arange(0.5*DR, (NR+0.5)*DR, DR) # bin centers
	if len(r_plot) > NR:
		if r_plot[-1] > r_max:
			r_plot = r_plot[:-1]
		else:
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

print 'Galaxies with suppression: ' + str(gal_ids[suppress_mask])
with h5py.File(results_dir+'suppressed_galaxies.h5', 'a') as f:
	f.create_dataset(model+'_'+snap, data=np.array(gal_ids[suppress_mask]))