import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
import h5py
import caesar
import pygad as pg
from yt import YTArray, YTQuantity
import sys
sys.path.append('/home/sapple/tools/')
from projection import recentre_pos_and_vel, compute_rotation_to_vec, rotate

def make_image(posx, posy, Npixels, mass, filename):
	xmin = -20.
	xmax = 20.
	im,xedges,yedges=np.histogram2d(posx,posy,bins=(Npixels,Npixels),weights=mass)
	im=im/((xmax-xmin)/float(Npixels))**2
	extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
	v_min = np.min(np.log10(im[im>0]))
	v_max = np.max(np.log10(im[im>0]))

	plt.imshow(np.log10(im.transpose()+0.0001),extent=extent, interpolation='nearest',cmap='magma',
				vmin=v_min,vmax=v_max, origin="lower")
	plt.savefig(filename)
	plt.clf()

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
	return np.sqrt(grav_energy - 0.5*(sigma_vel**2))

def center_of_quantity(quantity, weight):
	if len(quantity.shape) == 2:
		weight =  np.transpose(np.array([weight,]*len(quantity[0])))
	return np.sum(quantity*weight, axis=0) / np.sum(weight, axis=0)

def vrot_los(pos, vel, mass, vec, NR, DR, weight, r_plot, filename):
	axis, angle = compute_rotation_to_vec(pos[:, 0], pos[:, 1], pos[:, 2], vel[:, 0], vel[:, 1], vel[:, 2], mass, vec)
	pos_edgeon = rotate(pos[:, 0], pos[:, 1], pos[:, 2], axis, angle)
	vel_edgeon = rotate(vel[:, 0], vel[:, 1], vel[:, 2], axis, angle)
	r_edgeon = np.sqrt(pos_edgeon[0]**2 + pos_edgeon[1]**2)

	vrot = np.zeros(NR)
	for j in range(0,NR):
		edgeon_mask = (r_edgeon >= j*DR)*(r_edgeon < (j+1)*DR)
		if True in edgeon_mask:
			# find the rotational velocity in line of sight
			vrot[j] = center_of_quantity(vel_edgeon[2][edgeon_mask], weight[edgeon_mask])

	plt.plot(r_plot,vrot)
	plt.xlabel('R (kpc)')
	plt.ylabel('Vrot los (km/s)')
	plt.savefig(filename)
	plt.clf()

	return np.max(vrot)


if __name__ == '__main__':

	model = 'm50n512'
	wind = 's50j7k'
	snap = '151'
	results_dir = '/home/sapple/simba_sizes/sigma_rot/'
	data_dir = '/home/rad/data/'+model+'/'+wind+'/'

	s = pg.Snap(data_dir+'snap_'+model+'_'+snap+'.hdf5')
	sim =  caesar.load(data_dir+'Groups/'+model+'_'+snap+'.hdf5')

	h = sim.simulation.hubble_constant
	redshift = sim.simulation.redshift
	G = sim.simulation.G.in_units('km**3/(Msun*s**2)')

	ssfr_min = 0. # look at literature for this range
	sm_min = 4.6e9
	sfr_min = -1.5
	factor = 2. # do profile up to twice half mass radius
	NR = 10
	Npixels = 100

	edge_vec = np.array([0, 1, 0]) # for edge-on projection

	gal_sm = np.array([i.masses['stellar'].in_units('Msun') for i in sim.central_galaxies])
	gal_bm = np.array([i.masses['bh'].in_units('Msun') for i in sim.central_galaxies])
	gal_pos = np.array([i.pos.in_units('kpc') for i in sim.central_galaxies])
	gal_rad = np.array([i.radii['total_half_mass'].in_units('kpc') for i in sim.central_galaxies])
	gal_vel = np.array([i.vel.in_units('km/s') for i in sim.central_galaxies])
	gal_sfr = np.array([i.sfr.in_units('Msun/yr') for i in sim.central_galaxies])
	gal_ssfr = gal_sfr / gal_sm

	len_mask = np.array([len(i.glist) > 256 for i in sim.central_galaxies])
	sfr_mask = (np.log10(gal_sfr) > sfr_min) 
	sm_mask = (gal_sm > sm_min)
	gal_ids = np.arange(0, len(sim.central_galaxies))[sfr_mask*sm_mask*len_mask]

	#gas_v = readsnap(snapfile,'vel','gas',units=0,suppress=1)
	#gas_pos = readsnap(snapfile, 'pos', 'gas', suppress=1, units=1) / (h*(1.+redshift)) # in kpc
	#gas_mass = readsnap(snapfile, 'mass', 'gas', suppress=1, units=1) / h

	gas_v = s.gas['vel'].in_units_of('km/s')
	gas_pos = s.gas['pos'].in_units_of('ckpc') / (1+redshift) # in kpc
	gas_mass = s.gas['mass'].in_units_of('Msol')
	abundance_h1 = s.gas['NeutralHydrogenAbundance']
	abundance_h2 = s.gas['fh2']

	sigv_h1 = np.zeros(len(gal_ids)) # velocity dispersion centered on h1
	sigv_h2 = np.zeros(len(gal_ids)) # velocity dispersion centered on h2
	vrot_grav_h1 = np.zeros(len(gal_ids)) # using velocity dispersion from h1
	vrot_grav_h2 = np.zeros(len(gal_ids)) # using velocity dispersion from h2
	vrot_h1 = np.zeros(len(gal_ids))
	vrot_h2 = np.zeros(len(gal_ids))

	for i in range(len(gal_ids)):
		glist = sim.central_galaxies[gal_ids[i]].glist

		h1 = abundance_h1[glist]
		h2 = abundance_h2[glist]
		mass = gas_mass[glist]
		r_max = factor*gal_rad[gal_ids[i]]

		DR = r_max / NR
		r_plot = np.arange(0.5*DR, (NR+0.5)*DR, DR) # bin centers
		if len(r_plot) > NR:
			r_plot = r_plot[1:]

		# do this all for h1
		pos = gas_pos[glist]
		vel = gas_v[glist]
		pos -= center_of_quantity(pos, h1*mass)
		vel -= center_of_quantity(vel, h1*mass)
		mask = np.linalg.norm(pos, axis=1) <= r_max
		mass_within_r = YTQuantity(np.sum(mass[mask]), 'Msun')

		sigv_h1[i] = sigma_vel(vel[mask])
		# get mass within r_max or within position of maximum velocity?
		vrot_grav_h1[i] = vrot_gravity(mass_within_r, YTArray(r_max, 'kpc').in_units('km'), YTQuantity(sigv_h1[i], 'km/s'), G)

		plot_name = results_dir+'profiles/vrot_h1_profile_gal_'+str(i)+'.png'
		vrot_h1[i] = vrot_los(pos[mask], vel[mask], mass[mask], edge_vec, NR, DR, (h1*mass)[mask], r_plot, plot_name)


		# do this all for h2
		pos = gas_pos[glist]
		vel = gas_v[glist]
		pos -= center_of_quantity(pos, h2*mass)
		vel -= center_of_quantity(vel, h2*mass)
		mask = np.linalg.norm(pos, axis=1) <= r_max
		mass_within_r = YTQuantity(np.sum(mass[mask]), 'Msun')

		sigv_h2[i] = sigma_vel(vel[mask])
		# get mass within r_max or within position of maximum velocity?
		vrot_grav_h2[i] = vrot_gravity(mass_within_r, YTArray(r_max, 'kpc').in_units('km'), YTQuantity(sigv_h2[i], 'km/s'), G)

		plot_name = results_dir+'profiles/vrot_h2_profile_gal_'+str(i)+'.png'
		vrot_h2[i] = vrot_los(pos[mask], vel[mask], mass[mask], edge_vec, NR, DR, (h2*mass)[mask], r_plot, plot_name)

	with h5py.File(results_dir+'full_kinematic_profiles.h5', 'a') as f:
		f.create_dataset('sigmav_h1', data=np.array(sigv_h1))
		f.create_dataset('sigmav_h2', data=np.array(sigv_h2))
		f.create_dataset('vrot_gravity_h1', data=np.array(vrot_grav_h1))
		f.create_dataset('vrot_gravity_h2', data=np.array(vrot_grav_h2))
		f.create_dataset('vrot_h1', data=np.array(vrot_h1))
		f.create_dataset('vrot_h2', data=np.array(vrot_h2))
		f.create_dataset('gal_ids', data=np.array(gal_ids))

