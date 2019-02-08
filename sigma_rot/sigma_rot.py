import numpy as np
import caesar
from readgadget import readsnap

import sys
sys.path.append('/home/sapple/tools/')
from projection import recentre_pos_and_vel, compute_rotation_to_vec, rotate

def sigma_vel(velocities):
	"""
	Find the velocity dispersion of an array of velocities.
	"""
	mean_vel = np.mean(velocities, axis=0)
	dev = (velocities - mean_vel)**2
	sigma = np.sqrt(sum(dev)/len(dev))
	return (np.sqrt(sigma[0]**2 + sigma[1]**2 + sigma[2]**2))

def sigma_vel_los(velocities):
	"""
	Find the velocity dispersion along a line of sight
	"""
	mean_vel = np.mean(velocities)
	dev = (velocities - mean_vel)**2
	return np.sqrt(sum(dev)/len(dev))

def vrot_projected(velocities, positions, mass, r_max, vec):
	# 1. compute center of mass for particles within a given radius and shift all particles
	posx, posy, posz, vx_cm, vy_cm, vz_cm = recentre_pos_and_vel(positions[:, 0], positions[:, 1], positions[:, 2], 
																velocities[:, 0], velocities[:, 1], velocities[:, 2], mass, r_max)
	# 2. get axis and angle of rotation
	axis, angle = compute_rotation_to_vec(posx, posy, posz, vx_cm, vy_cm, vz_cm, mass, vec)
	# 3. rotate velocities to line of sight direction
	v_rotated = rotate(velocities[:, 0] - vx_cm, velocities[:, 1] - vy_cm, velocities[:, 2] - vz_cm, axis, angle) # these velocity components are projected onto the y axis (edge on)
	# now we have rotated gas particle positions for line of sight velocities
	# the line of sight component is the y direction (from vec)
	los = np.where(vec == 1)[0]
	return v_rotated[los]

def vrot_gravity(mass, r, sigma_vel, G):
	"""
	Find the rotational velocity from total energy
	"""
	grav_energy =  G*mass / r 
	return np.sqrt(grav_energy - 0.5*sigma_vel)


model = 'm50n512'
wind = 's50j7k'
snap = '151'

results_dir = '/home/sapple/simba_sizes/sigma_rot/'
data_dir = '/home/rad/data/'+model+'/'+wind+'/'
snapfile = data_dir+'snap_'+model+'_'+snap+'.hdf5'

sim =  caesar.load(data_dir+'Groups/'+model+'_'+snap+'.hdf5')

h = sim.simulation.hubble_constant
redshift = sim.simulation.redshift
G = sim.simulation.G.in_units('km**3/(Msun*s**2)')

sfr_min = 0.5 # look at literature for this range
mass_lower = 4.6e9
factor = 2.

DR = 1. # kpc

edge_vec = np.array([0, 1, 0]) # for edge-on projection
face_vec = np.array([0, 0, 1]) # for face-on projection

gal_sm = np.array([i.masses['stellar'].in_units('Msun') for i in sim.central_galaxies])
gal_bm = np.array([i.masses['bh'].in_units('Msun') for i in sim.central_galaxies])
gal_pos = np.array([i.pos.in_units('kpc') for i in sim.central_galaxies])
gal_rad = np.array([i.radii['total_half_mass'].in_units('kpc') for i in sim.central_galaxies])
gal_vel = np.array([i.vel.in_units('km/s') for i in sim.central_galaxies])
gal_sfr = np.array([i.sfr.in_units('Msun/yr') for i in sim.central_galaxies])

len_mask = np.array([len(i.glist) > 256 for i in sim.central_galaxies])
gal_ids = np.arange(0, len(sim.central_galaxies))[(gal_sfr > sfr_min)*len_mask]

gas_v = readsnap(snapfile,'vel','gas',units=0,suppress=1)
gas_pos = readsnap(snapfile, 'pos', 'gas', suppress=1, units=1) / (h*(1.+redshift)) # in kpc
gas_mass = readsnap(snapfile, 'mass', 'gas', suppress=1, units=1) / h

sigv_faceon = [] # velocity dispersion
vrot_edgeon = []
vrot_grav = []

for i in range(len(gal_ids)):
	glist = sim.central_galaxies[gal_ids[i]].glist

	pos = gas_pos[glist] - gal_pos[gal_ids[i]]
	vel = gas_v[glist]
	mass = gas_mass[glist]
	r_max = factor*gal_rad[gal_ids[i]]
	NR = int(r_max/DR)+1

	# 1. compute center of mass for particles within a given radius and shift all particles
	posx, posy, posz, vx_cm, vy_cm, vz_cm = recentre_pos_and_vel(pos[:, 0], pos[:, 1], pos[:, 2], 
																vel[:, 0], vel[:, 1], vel[:, 2], mass, r_max)

	# 2. get axis and angle of rotation
	axis, angle = compute_rotation_to_vec(posx, posy, posz, vx_cm, vy_cm, vz_cm, mass, face_vec)
	# 3. rotate velocities to line of sight direction
	pos_faceon = rotate(posx, posy, posz, axis, angle)
	vel_faceon = rotate(vx_cm, vy_cm, vz_cm, axis, angle)
	rfaceon = np.sqrt(pos_faceon[0]**2 + pos_faceon[1]**2)

	# 2. get axis and angle of rotation
	axis, angle = compute_rotation_to_vec(posx, posy, posz, vx_cm, vy_cm, vz_cm, mass, edge_vec)
	# 3. rotate velocities to line of sight direction
	pos_edgeon = rotate(posx, posy, posz, axis, angle)
	vel_edgeon = rotate(vx_cm, vy_cm, vz_cm, axis, angle)
	#redgeon = np.sqrt(pos_edgeon[0]**2 + pos_edgeon[1]**2)
	
	sigv_profile = np.zeros(NR)
	vrot_g_profile = np.zeros(NR)
	vrot_l_profile = np.zeros(NR)

	for j in range(0,NR):
		shell_mask = (rfaceon >= j*DR)*(rfaceon < (j+1)*DR)
		inner_mask = (rfaceon <= (j+1)*DR)

		# get the velocity dispersion within the face-on radial bin
		sigv_profile[j] = sigma_vel_los(vel_faceon[shell_mask][2])

		# find the rotational velocity from gravitational energy
		mass_within_r = np.sum(mass[inner_mask])
		r_km = (j+1)*DR /(3.086e13*1.e3)
		vrot_g_profile[j] = vrot_gravity(mass_within_r, r_km, sigv3d[i][j], G)
		
		# find the rotational velocity in line of sight
		vrot_l_profile[j] = np.max(np.abs(vel_edgeon[1][shell_mask]))

	sigv_faceon.append(sigv_faceon)
	vrot_grav.append(vrot_g_profile)
	vrot_edgeon.append(vrot_l_profile)
	
# need to check all these projections and rotations etc.