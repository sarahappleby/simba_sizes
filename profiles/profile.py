import numpy as np 
import matplotlib.pyplot as plt
import h5py
import sys
import os

sys.path.append('/home/sapple/tools/')
from projection import *

from readgadget import readsnap
import caesar
from astropy.cosmology import FlatLambdaCDM

from profile_methods import *

ssfr_min = 3.e-10
sm_min = 5.e9
age_min = 50.
age_max = 150.
mass_loss = 1.18
time = 100.e6
do_ages = True
bh_center = True

# for making profiles
n = 10
factor = 2.
dr = factor / n
vec = np.array([0, 0, 1]) # face-on projection to collapse the z direction
bins = np.arange(0., factor, dr)
area = np.array([np.pi*(i+dr)**2 -   np.pi* i**2 for i in bins])

model = sys.argv[1]
wind = sys.argv[2]
snap = sys.argv[3]

if bh_center:
	results_dir = '/home/sapple/simba_sizes/profiles/'+model+'/snap_'+snap+'/bh_centered/'
else:
	results_dir = '/home/sapple/simba_sizes/profiles/'+model+'/snap_'+snap+'/scom_centered/'
if not os.path.exists(results_dir):
	os.makedirs(results_dir)
	os.makedirs(results_dir + 'images/')
	os.makedirs(results_dir + 'profiles/')

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

use_star_m = np.array([np.array([]) for i in range(n)])
use_young_star = np.array([np.array([]) for i in range(n)])
use_gas_sfr = np.array([np.array([]) for i in range(n)])
use_gas_h1 = np.array([np.array([]) for i in range(n)])
use_gas_h2 = np.array([np.array([]) for i in range(n)])

no_gals_star = np.zeros(n)
no_gals_gas = np.zeros(n)

for i in range(len(gal_ids)):

	print 'Galaxy ' +str(gal_ids[i])
	slist = sim.galaxies[gal_ids[i]].slist
	glist = sim.galaxies[gal_ids[i]].glist
	r_max = factor*gal_rad[gal_ids[i]]

	"""
	Find stellar ages:
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

	else:
		ages_mask = np.array([True for j in range(len(slist))])
	"""
	Get star particle data and correct for bh center or star com
	"""
	pos = star_pos[slist]
	vel = star_vels[slist]
	mass = star_mass[slist]
	if not bh_center:
		pos -= center_of_quantity(pos, mass)
	else:
		pos -= bh_pos[sim.galaxies[gal_ids[i]].bhlist[0]]
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
	Rotate to disk plane
	"""
	axis, angle = compute_rotation_to_vec(pos[:, 0], pos[:, 1], pos[:, 2], vel[:, 0], vel[:, 1], vel[:, 2], mass, vec)
	pos[:, 0], pos[:, 1], pos[:, 2] = rotate(pos[:, 0], pos[:, 1], pos[:, 2], axis, angle)
	plot_name = results_dir + 'images/gal_'+str(gal_ids[i]) + '_stars.png'
	make_image(pos[:, 0], pos[:, 1], np.log10(mass), plot_name, r_max)
	r = np.linalg.norm(pos[:, [0, 1]], axis=1) / r_max*factor
	"""
	Bin by radial distance
	"""
	digitize = np.digitize(r, bins) - np.ones(len(r), dtype='int64')
	binned_m = [mass[digitize == j] for j in range(0, len(bins))]
	binned_young_star = [ages_mask[digitize == j] for j in range(0, len(bins))]
	"""
	Add into main halo array
	"""
	use_star_m = np.array([np.append(use_star_m[i], binned_m[i]) for i in range(n)])
	use_young_star = np.array([np.append(use_young_star[i], binned_young_star[i]) for i in range(n)])
	no_gals_star += np.array([1. if len(binned_m[j]) > 0. else 0. for j in range(n)])



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
	"""
	Rotate to disk plane
	"""
	axis, angle = compute_rotation_to_vec(pos[:, 0], pos[:, 1], pos[:, 2], vel[:, 0], vel[:, 1], vel[:, 2], mass, vec)
	pos[:, 0], pos[:, 1], pos[:, 2] = rotate(pos[:, 0], pos[:, 1], pos[:, 2], axis, angle)
	plot_name = results_dir + 'images/gal_'+str(gal_ids[i]) + '_gas.png'
	make_image(pos[:, 0], pos[:, 1], np.log10(mass), plot_name, r_max)
	r = np.linalg.norm(pos[:, [0, 1]], axis=1) / r_max
	"""
	Bin by radial distance
	"""
	digitize = np.digitize(r, bins) - np.ones(len(r), dtype='int64')
	binned_sfr = [sfr[digitize == j] for j in range(0, len(bins))]
	binned_h1 = [h1[digitize == j] for j in range(0, len(bins))]
	binned_h2 = [h2[digitize == j] for j in range(0, len(bins))]
	"""
	Add into main halo array
	"""
	use_gas_sfr = np.array([np.append(use_gas_sfr[i], binned_sfr[i]) for i in range(n)])
	use_gas_h1 = np.array([np.append(use_gas_h1[i], binned_h1[i]) for i in range(n)])
	use_gas_h2 = np.array([np.append(use_gas_h2[i], binned_h2[i]) for i in range(n)])
	no_gals_gas += np.array([1. if len(binned_sfr[j]) > 0. else 0. for j in range(n)])



sm_profile = np.array([np.sum(j) for j in  use_star_m]) / (no_gals_star*area)
star_sfr_profile = (np.array([np.sum(j) for j in use_young_star * use_star_m]) * mass_loss / time) / (no_gals_star*area)
star_ssfr_profile = sfr_profile / sm_profile
gas_sfr_profile = np.array([np.sum(j) for j in use_gas_sfr]) / (no_gals_gas*area)
gas_sfr_h1_profile = np.array([np.sum(j) for j in use_gas_h1]) / (no_gals_gas*area)
gas_sfr_h2_profile = np.array([np.sum(j) for j in use_gas_h2]) / (no_gals_gas*area)
gas_ssfr_profile = gas_sfr_profile / sm_profile

with h5py.File(results_dir+'profiles.h5', 'a') as f:
	f.create_dataset('gas_sfr_profile', data=np.array(use_gas_sfr))
	f.create_dataset('h1_profile', data=np.array(use_gas_h1))
	f.create_dataset('h2_profile', data=np.array(use_gas_h2))
	f.create_dataset('sm_profile', data=np.array(use_star_m))
	f.create_dataset('young_stars', data=np.array(use_young_star))
	f.create_dataset('gal_ids', data=np.array(gal_ids))
	f.create_dataset('no_gals_star', data=np.array(no_gals_star))
	f.create_dataset('no_gals_gas', data=np.array(no_gals_gas))


plt.semilogy(bins+(dr*0.5), sm_profile, marker='.', markersize=2)
plt.xlabel('R half *')
plt.ylabel('M* surface density (Msun/R half bin^2)')
plt.xlim(0, )
plt.title(str(len(gal_ids))+' galaxies')
plt.savefig(results_dir+'sm_profile.png')
plt.clf()

plt.plot(bins+(dr*0.5), star_sfr_profile, marker='.', markersize=2)
plt.xlabel('R half *')
plt.ylabel('SFR surface density (Msun/yr /R half bin^2)')
plt.xlim(0, )
plt.fill_between(bins+0.1, star_sfr_lower, star_sfr_higher, facecolor='blue', alpha=0.3)
plt.title(str(len(gal_ids))+' galaxies')
plt.savefig(results_dir+'star_sfr_profile.png')
plt.clf()

plt.plot(bins+(dr*0.5), np.log10(star_ssfr_profile), marker='.', markersize=2)
plt.xlabel('R half *')
plt.ylabel('log sSFR (yr^-1)')
plt.xlim(0, )
plt.title(str(len(gal_ids))+' galaxies')
plt.savefig(results_dir+'star_ssfr_profile.png')
plt.clf()

plt.plot(bins+(dr*0.5), gas_sfr_profile, marker='.', markersize=2)
plt.xlabel('R half *')
plt.ylabel('SFR surface density (Msun/yr /R half bin^2)')
plt.xlim(0, )
plt.title(str(len(gal_ids_sm))+' galaxies')
plt.savefig(results_dir+'gas_sfr_profile.png')
plt.clf()

plt.plot(bins+(dr*0.5), np.log10(gas_ssfr_profile), marker='.', markersize=2)
plt.xlabel('R half *')
plt.ylabel('log sSFR (yr^-1)')
plt.xlim(0, )
plt.title(str(len(gal_ids_sm))+' galaxies')
plt.savefig(results_dir+'gas_ssfr_profile.png')
plt.clf()

plt.plot(bins+(dr*0.5), h1_profile, marker='.', markersize=2)
plt.xlabel('R half *')
plt.ylabel('HI surface density (Msun /R half bin^2)')
plt.xlim(0, )
plt.title(str(len(gal_ids_sm))+' galaxies')
plt.savefig(results_dir+'h1_profile.png')
plt.clf()

plt.plot(bins+(dr*0.5), h2_profile, marker='.', markersize=2)
plt.xlabel('R half *')
plt.ylabel('HII surface density (Msun /R half bin^2)')
plt.xlim(0, )
plt.title(str(len(gal_ids_sm))+' galaxies')
plt.savefig(results_dir+'h2_profile.png')
plt.clf()

