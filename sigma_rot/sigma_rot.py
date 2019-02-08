# need to get v_circ

import numpy as np
import caesar
from readgadget import readsnap

import sys
sys.path.append('/home/sapple/tools/')
from projection import *

model = 'm50n512'
wind = 's50j7k'
snap = '151'

results_dir = '/home/sapple/simba_sizes/sigma_rot/'
data_dir = '/home/rad/data/'+model+'/'+wind+'/'
snapfile = data_dir+'snap_'+model+'_'+snap+'.hdf5'

sim =  caesar.load(data_dir+'Groups/'+model+'_'+snap+'.hdf5')

h = sim.simulation.hubble_constant
redshift = sim.simulation.redshift

sfr_min = 0.5 # look at literature for this range
mass_lower = 4.6e9
factor = 5.

vec = np.array([0, 1, 0]) # for edge-on projection

gal_sm = np.array([i.masses['stellar'].in_units('Msun') for i in sim.central_galaxies])
gal_bm = np.array([i.masses['bh'].in_units('Msun') for i in sim.central_galaxies])
gal_pos = np.array([i.pos.in_units('kpc') for i in sim.central_galaxies])
gal_rad = np.array([i.radii['total_half_mass'].in_units('kpc') for i in sim.central_galaxies])
gal_vel = np.array([i.vel.in_units('km/s') for i in sim.central_galaxies])

gal_sfr = np.array([i.sfr.in_units('Msun/yr') for i in sim.central_galaxies])
gal_ids = np.arange(0, len(sim.central_galaxies))[gal_sfr > sfr_min]

gas_v = readsnap(snapfile,'vel','gas',units=0,suppress=1)
gas_pos = readsnap(snapfile, 'pos', 'gas', suppress=1, units=1) / (h*(1.+redshift)) # in kpc
gas_vels = readsnap(snapfile, 'vel', 'gas', suppress=1, units=0)
gas_mass = readsnap(snapfile, 'mass', 'gas', suppress=1, units=1) / (h*10**7.)

sigv3d = np.zeros(len(gal_ids)) # velocity dispersion
vrot = np.zeros(len(gal_ids))

j = 0
for i in gal_ids:
	glist = sim.central_galaxies[i].glist

	if len(glist) < 256:
		print 'Galaxy ' + str(i) + ': not enough particles'
		continue

	# get the velocity dispersion
	svgal = np.array([gas_v[k] for k in glist])
	svcent = np.mean(svgal,axis=0)
	svgal = (svgal-svcent)**2
	sigv = np.sqrt(sum(svgal)/len(svgal))
	sigv3d[j] = (np.sqrt(sigv[0]**2 + sigv[1]**2 + sigv[2]**2)) # done

	# kassin et al 2014 - rotational velocity from edge-on line of sight velocities

	# do projection for rotational velocity
	x = gas_pos[glist][:, 0] - gal_pos[i][0]
	y = gas_pos[glist][:, 1] - gal_pos[i][1]
	z = gas_pos[glist][:, 2] - gal_pos[i][2]
	vx = gas_vels[glist][:, 0]
	vy = gas_vels[glist][:, 1]
	vz = gas_vels[glist][:, 2]
	mass = gas_mass[glist]

	r_max = factor*gal_rad[i]
	# 1. compute center of mass for particles within a given radius and shift all particles
	posx, posy, posz, vx_cm, vy_cm, vz_cm = recentre_pos_and_vel(x, y, z, vx, vy, vz, mass, r_max)
	# 2. get axis and angle of rotation
	axis, angle = compute_rotation_to_vec(posx, posy, posz, vx_cm, vy_cm, vz_cm, mass, vec)
	# 3. rotate velocities to line of sight direction
	rel_vx, rel_vy, rel_vz = rotate(vx - vx_cm, vy - vy_cm, vz - vz_cm, axis, angle) # these velocity components are projected onto the y axis (edge on)
	# now we have rotated gas particle positions for line of sight velocities
	# the line of sight component is the y direction (from vec)
	vrot[j] = np.max(np.abs(rel_vy))

	j +=1