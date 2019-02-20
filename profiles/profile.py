import sys
sys.path.append('/home/sapple/tools/')
from projection import *

import numpy as np 
from readgadget import readsnap
import caesar
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
import h5py

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
snap = '078'

results_dir = '/home/sapple/simba_sizes/profiles/snap_'+snap+'/'
data_dir = '/home/rad/data/'+model+'/'+wind+'/'
snapfile = data_dir+'snap_'+model+'_'+snap+'.hdf5'

sim =  caesar.load(data_dir+'Groups/'+model+'_'+snap+'.hdf5')

h = sim.simulation.hubble_constant
redshift = sim.simulation.redshift
cosmo = FlatLambdaCDM(H0=100*h, Om0=sim.simulation.omega_matter, Ob0=sim.simulation.omega_baryon, Tcmb0=2.73)
thubble = cosmo.age(redshift).value # in Gyr

# for sample of galaxies and stars
ssfr_min = 1.e-10
mass_lower = 4.6e9
age_min = 50.
age_max = 100.

# for making profiles
NR = 10
factor = 2.
vec = np.array([0, 0, 1]) # face-on projection to collapse the z direction

# for plotting
xmin = -20.
xmax = 20.
Npixels = 100

# load in the caesar galaxy data to make an initial cut of star forming galaxies
gal_sm = np.array([i.masses['stellar'].in_units('Msun') for i in sim.central_galaxies])
gal_bm = np.array([i.masses['bh'].in_units('Msun') for i in sim.central_galaxies])
gal_pos = np.array([i.pos.in_units('kpc') for i in sim.central_galaxies])
gal_rad = np.array([i.radii['total_half_mass'].in_units('kpc') for i in sim.central_galaxies])
gal_sfr = np.array([i.sfr.in_units('Msun/yr') for i in sim.central_galaxies])
gal_ssfr = gal_sfr / gal_sm

sm_mask = gal_sm > sm_min
ssfr_mask = gal_ssfr > ssfr_min
gal_ids = np.arange(0, len(sim.central_galaxies))[ssfr_mask*sm_mask]

# load in star particle data with readsnap
star_pos = readsnap(snapfile, 'pos', 'star', suppress=1, units=1) / (h*(1.+redshift)) # in kpc
star_mass = readsnap(snapfile, 'mass', 'star', suppress=1, units=1) / h
star_vels = readsnap(snapfile, 'vel', 'star', suppress=1, units=0)
star_tform = readsnap(snapfile, 'age', 'star', suppress=1, units=1) # expansion times at time of formation

ssfr_profiles = np.zeros((len(gal_ids), NR))
sm_profiles = np.zeros((len(gal_ids), NR))

len_mask = np.ones(len(gal_ids))

for i in range(len(gal_ids)):
	slist = sim.central_galaxies[gal_ids[i]].slist
	
	# find ages of galaxy stars
	tform = star_tform[slist]
	ages = tage(cosmo, thubble, tform) * 1000. # get star ages in Myr
	ages_mask = (ages > age_min) & (ages < age_max)

	# filter for galaxies with more than 256 new star particles for the profiles
	if len(ages[]ages_mask) < 256:
		len_mask[i] = 0.
		continue

	pos = gas_pos[slist]
	vel = gas_v[slist]
	mass = star_mass[slist]
	r_max = factor*gal_rad[gal_ids[i]]

	DR = r_max / NR
	r_plot = np.arange(0.5*DR, (NR+0.5)*DR, DR) # bin centers
	if len(r_plot) > NR:
		r_plot = r_plot[1:]

	pos_ssfr = center_of_quantity(pos, ssfr)

	posx, posy, posz, vx, vy, vz = recentre_pos_and_vel(x, y, z, vx, vy, vz, mass, r_max)
			
	# second center of mass recentering
	filter_rad = (np.sqrt(posx**2+posy**2+posz**2) < r_max)
	posx, posy, posz, vx, vy, vz = recentre_pos_and_vel(posx[filter_rad], posy[filter_rad], posz[filter_rad], 
														vx[filter_rad], vy[filter_rad], vz[filter_rad], mass[filter_rad], r_max)
	mass = mass[filter_rad]
	ages_mask = ages_mask[filter_rad]

	# 2. get axis and angle of rotation
	axis, angle = compute_rotation_to_vec(posx, posy, posz, vx, vy, vz, mass, vec)
	# 3. rotate positions and velocities
	posx, posy, posz = rotate(posx, posy, posz, axis, angle)

	# make image of face-on galaxy
	filter_rad=(posx>xmin)*(posx<xmax)*(posy>xmin)*(posy<xmax)
	im,xedges,yedges=np.histogram2d(posx[filter_rad],posy[filter_rad],bins=(Npixels,Npixels),weights=mass[filter_rad])
	im=im/((xmax-xmin)/float(Npixels))**2 #gives surface density
	extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

	v_min = np.min(np.log10(im[im>0]))
	v_max = np.max(np.log10(im[im>0]))

	plt.imshow(np.log10(im.transpose()+0.0001),extent=extent, interpolation='nearest',cmap='magma',vmin=v_min,vmax=v_max, origin="lower")
	plt.title('Mass: ' + str.format("{0:.6g}", float(gal_sm[i])) + ' Msun')
	plt.savefig(results_dir + 'image_stars_gal_'+str(i)+ '.png')
	plt.clf()

	# make profile of total mass of stars with ages 50Myr - 100Myr
	mass_new = mass[ages_mask]*10**7
	r = np.sqrt(posx*posx + posy*posy)[ages_mask]

	for j in range(0,NR):
		mask = (r >= j*DR)*(r < (j+1)*DR)
		sm_profiles[i][j] = np.sum(mass_new[mask])

	# divide by 50Myr to get SFR in last 50Myr with 18% instantaneous mass loss
	ssfr_profiles[i] = sm_profiles[i]*1.18 / 50.e6

	# plotting profiles
	plt.semilogy(r_plot,sm_profiles[i])
	plt.xlabel('R (kpc)')
	plt.ylabel('M_* (Msun)')
	plt.savefig(results_dir+'profiles/sm_profile_gal_'+str(i)+'.png')
	plt.clf()

	plt.plot(r_plot,ssfr_profiles[i])
	plt.xlabel('R (kpc)')
	plt.ylabel('sSFR (Msun/yr)')
	plt.savefig(results_dir+'profiles/ssfr_profile_gal_'+str(i)+'.png')
	plt.clf()

with h5py.File(results_dir+model+'_'+snap+'_profiles.h5', 'a') as f:
	f.create_dataset('ssfr_profiles', data=np.array(ssfr_profiles))
	f.create_dataset('sm_profiles', data=np.array(sm_profiles))
	f.create_dataset('gal_ids', data=np.array(gal_ids))
	f.create_dataset('gal_sm', data=np.array(gal_sm))
	f.create_dataset('gal_bm', data=np.array(gal_bm))
	f.create_dataset('gal_sfr', data=np.array(gal_sfr))
f.close()
