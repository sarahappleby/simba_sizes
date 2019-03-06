import sys
sys.path.append('/home/sapple/tools/')
from projection import *

import numpy as np 
from readgadget import readsnap
import caesar
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
import h5py

def make_image(posx, posy, weight, filename, Npixels=100):
	xmin = -20.
	xmax = 20.
	im,xedges,yedges=np.histogram2d(posx,posy,bins=(Npixels,Npixels),weights=weight)
	im=im/((xmax-xmin)/float(Npixels))**2
	extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
	v_min = np.min(np.log10(im[im>0]))
	v_max = np.max(np.log10(im[im>0]))

	plt.imshow(np.log10(im.transpose()+0.0001),extent=extent, interpolation='nearest',cmap='magma',
				vmin=v_min,vmax=v_max, origin="lower")
	plt.savefig(filename)
	plt.clf()

def center_of_quantity(quantity, weight):
	if len(quantity.shape) == 2:
		weight =  np.transpose(np.array([weight,]*len(quantity[0])))
	return np.sum(quantity*weight, axis=0) / np.sum(weight, axis=0)

def tage(cosmo,thubble,a):
	"""
	Find age of stars in Gyr from expansion time at time of formation
	"""
	return thubble-cosmo.age(1./a-1).value

model = 'm50n512'
wind = 's50j7k'
snap = '151'

results_dir = '/home/sapple/simba_sizes/profiles/snap_'+snap+'/'
data_dir = '/home/rad/data/'+model+'/'+wind+'/'
snapfile = data_dir+'snap_'+model+'_'+snap+'.hdf5'

# for sample of galaxies and stars
ssfr_min = 3.e-10
sm_min = 5.e9
age_min = 50.
age_max = 100.
mass_loss = 1.18
time = 50.e6
do_ages = True

# for making profiles
DR = 1.
factor = 2.
vec = np.array([0, 0, 1]) # face-on projection to collapse the z direction
bins = np.arange(0., 2., 0.2)

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

with h5py.File('/home/sapple/simba_sizes/sizes/data/'+model+'_'+wind+'_'+snap+'_data.h5', 'r') as f:
	gal_rad = f['halfmass'].value

# load in star particle data with readsnap
star_pos = readsnap(snapfile, 'pos', 'star', suppress=1, units=1) / (h*(1.+redshift)) # in kpc
star_mass = readsnap(snapfile, 'mass', 'star', suppress=1, units=1) / h
star_vels = readsnap(snapfile, 'vel', 'star', suppress=1, units=0)
star_tform = readsnap(snapfile, 'age', 'star', suppress=1, units=1) # expansion times at time of formation


ssfr_profiles = np.zeros((len(gal_ids), 10))
sm_profiles = np.zeros((len(gal_ids), 10))
sm_frac_profiles = np.zeros((len(gal_ids), 10))

len_mask = np.array([True for i in range(len(gal_ids))])

for i in range(len(gal_ids)):

	print 'Galaxy ' +str(gal_ids[i])
	slist = sim.galaxies[gal_ids[i]].slist
	if do_ages:
		print 'Find stellar ages'
		# find ages of galaxy stars
		tform = star_tform[slist]
		ages = tage(cosmo, thubble, tform) * 1000. # get star ages in Myr
		ages_mask = (ages > age_min) & (ages < age_max)

		# filter for galaxies with more than 256 new star particles for the profiles
		if len(ages[ages_mask]) < 256:
			len_mask[i] = False
			continue
	else:
		ages_mask = [True for j in range(len(slist))]
		
	pos = star_pos[slist]
	vel = star_vels[slist]
	mass = star_mass[slist]
	total_mass = np.sum(mass)
	r_max = factor*gal_rad[gal_ids[i]]

	# for the profile
	NR = int(round(r_max / DR))
	if NR < 5:
		len_mask[i] = False
		print 'Galaxy ' + str(gal_ids[i]) + ' too small (below 5 kpc)'
		continue

	pos -= center_of_quantity(pos, mass)
	vel -= center_of_quantity(vel, mass)
	r = np.linalg.norm(pos, axis=1)
	mass = mass[r < r_max ]
	vel = vel[r < r_max ]
	pos = pos[r < r_max ]
	ages_mask = np.array(ages_mask)[r < r_max]
	r = r[r < r_max]

	axis, angle = compute_rotation_to_vec(pos[:, 0], pos[:, 1], pos[:, 2], vel[:, 0], vel[:, 1], vel[:, 2], mass, vec)
	pos[:, 0], pos[:, 1], pos[:, 2] = rotate(pos[:, 0], pos[:, 1], pos[:, 2], axis, angle)
	plot_name = results_dir + 'images/gal_'+str(gal_ids[i]) + '.png'
	make_image(pos[:, 0], pos[:, 1], mass, plot_name)

	profile = np.zeros(NR)
	ages_profile = np.zeros(NR)
	# make profile of total mass of stars with ages 50Myr - 100Myr
	for j in range(0, NR):
		mask = (r >= j*DR)*(r < (j+1)*DR)
		profile[j] = np.sum(mass[mask])
		ages_profile[j] = np.sum(mass[mask*ages_mask])
		if (j==0):
			profile[j] /= np.pi*DR*DR
			ages_profile /= np.pi*DR*DR
		else:
			profile[j] /= np.pi*(DR*DR*(j+1)*(j+1) - DR*DR*j*j)
			ages_profile[j] /= np.pi*(DR*DR*(j+1)*(j+1) - DR*DR*j*j)
	# divide by 50Myr to get SFR in last 50Myr with 18% instantaneous mass loss
	ssfr_profile = ages_profile*mass_loss / time
	sm_frac = profile / total_mass

	r_plot = np.arange(0.5*DR, (NR+0.5)*DR, DR) # bin centers
	if len(r_plot) > NR:
		r_plot = r_plot[1:]

	plt.semilogy(r_plot,profile)
	plt.xlabel('R (kpc)')
	plt.ylabel('M* surface density (Msun/kpc^2)')
	plt.savefig(results_dir+'profiles/sm_profile_gal_'+str(gal_ids[i])+'.png')
	plt.clf()

	plt.plot(r_plot,ssfr_profile)
	plt.xlabel('R (kpc)')
	plt.ylabel('sSFR surface density (Msun/yr kpc^2)')
	plt.savefig(results_dir+'profiles/ssfr_profile_gal_'+str(gal_ids[i])+'.png')
	plt.clf()

	x = np.arange(0, NR) / r_max*factor
	digitized = np.digitize(x, bins)
	sm_bin = [profile[digitized == j] for j in range(1, len(bins) +1)]
	ssfr_bin = [ssfr_profile[digitized == j] for j in range(1, len(bins) +1)]
	sm_frac_bin = [sm_frac[digitized == j] for j in range(1, len(bins) +1)]
	
	ssfr_profiles[i] = np.array([np.sum(j) for j in ssfr_bin])
	sm_profiles[i] = np.array([np.sum(j) for j in sm_bin])
	sm_frac_profiles[i] = np.array([np.sum(j) for j in sm_frac_bin])

sm_profiles = sm_profiles[len_mask]
sm_frac_profiles = sm_frac_profiles[len_mask]
ssfr_profiles = ssfr_profiles[len_mask]
gal_ids = gal_ids[len_mask]

ssfr_lower = np.percentile(ssfr_profiles, 25, axis=0)
ssfr_higher = np.percentile(ssfr_profiles, 75, axis=0)
ssfr_median = np.percentile(ssfr_profiles, 50, axis=0)

sm_lower = np.percentile(sm_profiles, 25, axis=0)
sm_higher = np.percentile(sm_profiles, 75, axis=0)
sm_median = np.percentile(sm_profiles, 50, axis=0)

sm_frac_lower = np.percentile(sm_frac_profiles, 25, axis=0)
sm_frac_higher = np.percentile(sm_frac_profiles, 75, axis=0)
sm_frac_median = np.percentile(sm_frac_profiles, 50, axis=0)

plt.plot(bins+0.1, sm_median, marker='.', markersize=2)
plt.xlabel('R half *')
plt.ylabel('M* surface density (Msun/kpc^2)')
plt.fill_between(bins+0.1, sm_lower, sm_higher, facecolor='blue', alpha=0.3)
plt.savefig(results_dir+'sm_profile.png')
plt.clf()

plt.plot(bins+0.1, sm_frac_median, marker='.', markersize=2)
plt.xlabel('R half *')
plt.ylabel('M* surface density fraction (kpc^-2)')
plt.fill_between(bins+0.1, sm_frac_lower, sm_frac_higher, facecolor='blue', alpha=0.3)
plt.savefig(results_dir+'sm_frac_profile.png')
plt.clf()

plt.plot(bins+0.1, ssfr_median, marker='.', markersize=2)
plt.xlabel('R half *')
plt.ylabel('sSFR surface density (Msun /yr kpc^2)')
plt.fill_between(bins+0.1, ssfr_lower, ssfr_higher, facecolor='blue', alpha=0.3)
plt.savefig(results_dir+'ssfr_profile.png')
plt.clf()

with h5py.File(results_dir+model+'_'+snap+'_profiles.h5', 'a') as f:
	f.create_dataset('ssfr_profile', data=np.array(ssfr_profiles))
	f.create_dataset('sm_profile', data=np.array(sm_profiles))
	f.create_dataset('gal_ids', data=np.array(gal_ids))
	f.create_dataset('sm_frac_profiles', data=np.array(sm_frac_profiles))
