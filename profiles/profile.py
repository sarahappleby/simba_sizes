import sys
sys.path.append('/home/sapple/tools/')
from projection import *

import numpy as np 
from readgadget import readsnap
import caesar
from astropy.cosmology import FlatLambdaCDM
import matplotlib.pyplot as plt
import h5py

def tage(cosmo,thubble,a):
	"""
	Find age of stars in Gyr from expansion time at time of formation
	"""
	return thubble-cosmo.age(1./a-1).value

model = 'm50n512'
wind = 's50j7k'
snap = '078'

results_dir = '/home/sapple/simba_sizes/profiles/profiles/'
data_dir = '/home/rad/data/'+model+'/'+wind+'/'
snapfile = data_dir+'snap_'+model+'_'+snap+'.hdf5'

sim =  caesar.load(data_dir+'Groups/'+model+'_'+snap+'.hdf5')

h = sim.simulation.hubble_constant
z = sim.simulation.redshift
cosmo = FlatLambdaCDM(H0=100*h, Om0=sim.simulation.omega_matter, Ob0=sim.simulation.omega_baryon, Tcmb0=2.73)
thubble = cosmo.age(z).value # in Gyr

# for sample of galaxies and stars
sfr_min = 0.5 # look at literature for this range
mass_lower = 4.6e9
age_min = 50.
age_max = 100.

# for making profiles
disk_profile = True
spherical_profile = False
r_rot = 500.
DR = 1. # kpc
NR = 30
r_plot = np.arange(0,NR*DR,DR)

# for plotting
xmin = -20.
xmax = 20.
Npixels = 100

# load in the caesar galaxy data to make an initial cut of star forming galaxies
gal_sm = np.array([i.masses['stellar'].in_units('Msun') for i in sim.galaxies])
gal_bm = np.array([i.masses['bh'].in_units('Msun') for i in sim.galaxies])
gal_pos = np.array([i.pos.in_units('kpc') for i in sim.galaxies])
gal_rad = np.array([i.radius.in_units('kpc') for i in sim.galaxies])

gal_sfr = np.array([i.sfr.in_units('Msun/yr') for i in sim.galaxies])
gal_ids = np.arange(0, len(sim.galaxies))[gal_sfr > sfr_min]

# load in star particle data with readsnap
star_pos = readsnap(snapfile, 'pos', 'star', suppress=1, units=1) / (h*(1.+z)) # in kpc
star_mass = readsnap(snapfile, 'mass', 'star', suppress=1, units=1) / (h*10**7.)
star_vels = readsnap(snapfile, 'vel', 'star', suppress=1, units=0)
star_tform = readsnap(snapfile, 'age', 'star', suppress=1, units=1) # expansion times at time of formation

ssfr_profiles = np.zeros((len(gal_ids), NR))
sm_profiles = np.zeros((len(gal_ids), NR))

for i in gal_ids:
	slist = sim.galaxies[i].slist

	x = star_pos[slist][:, 0] - gal_pos[i][0]
	y = star_pos[slist][:, 1] - gal_pos[i][1]
	z = star_pos[slist][:, 2] - gal_pos[i][2]
	vx = star_vels[slist][:, 0]
	vy = star_vels[slist][:, 1]
	vz = star_vels[slist][:, 2]
	mass = star_mass[slist]

	# find ages of galaxy stars
	tform = star_tform[slist]
	ages = tage(cosmo, thubble, tform) * 1000. # get star ages in Myr
	ages_mask = (ages > age_min) & (ages < age_max)

	# filter for galaxies with more than 256 new star particles for the profiles
	if len(np.where(ages_mask == True)[0]) < 256:
		continue

	if spherical_profile:
		pass

	if disk_profile:
		vec = np.array([0, 0, 1]) # face-on projection to collapse the z direction

		# 1. compute center of mass for particles within a given radius and shift all particles
		r_max = 3.*gal_rad[i]
		posx, posy, posz, vx, vy, vz = recentre_pos_and_vel(x, y, z, vx, vy, vz, mass, r_max)
		filter_rad = (np.sqrt(posx**2+posy**2+posz**2) < 4.*gal_rad[i])
		# 2. get axis and angle of rotation
		axis, angle = compute_rotation_to_vec(posx[filter_rad], posy[filter_rad], posz[filter_rad],
											  vx[filter_rad], vy[filter_rad], vz[filter_rad], mass[filter_rad], vec)
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
		plt.savefig(results_dir + 'gal_image_'+str(i)+ '.png')
		plt.clf()

		# make profile of total mass of stars with ages 50Myr - 100Myr
		mass_new = mass[ages_mask]*10**7
		r = np.sqrt(posx*posx + posy*posy)[ages_mask]

		for j in range(0,NR):
			mask = (r >= j*DR)*(r < (j+1)*DR)
			sm_profiles[i][j] = np.sum(mass_new[mask])

		# divide by 50Myr to get SFR in last 50Myr
		ssfr_profiles[i] = sm_profiles[i] / 50.

		# plotting profiles
		plt.semilogy(r_plot,sm_profiles[i])
		plt.xlabel('R (kpc)')
		plt.ylabel('M_*')
		plt.savefig(results_dir+'gal_'+str(i)+'_sm_profile.png')
		plt.clf()

		plt.plot(r_plot,ssfr_profiles[i])
		plt.xlabel('R (kpc)')
		plt.ylabel('sSFR')
		plt.savefig(results_dir+'gal_'+str(i)+'_ssfr_profile.png')
		plt.clf()

with h5py.File(results_dir+model+'_'+snap+'_profiles.h5', 'a') as f:
	f.create_dataset('ssfr_profiles', data=np.array(ssfr_profiles))
	f.create_dataset('sm_profiles', data=np.array(sm_profiles))
	f.create_dataset('gal_ids', data=np.array(gal_ids))
	f.create_dataset('gal_sm', data=np.array(gal_sm))
	f.create_dataset('gal_bm', data=np.array(gal_bm))
	f.create_dataset('gal_sfr', data=np.array(gal_sfr))
f.close()
