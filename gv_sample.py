import caesar
import numpy as np

def sfms_line(x0, y0):
	a = 0.73
	b = -7.7
	y = x0*a + b
	if y < y0:
		return True
	else:
		return False

# in belfiore 18, fig 2:
sfr_min = -1.5
sm_min = 10.
sm_max = 11.5

model = 'm50n512'
wind = 's50j7k'
snap = '151'

results_dir = '/home/sapple/simba_sizes/sigma_rot/'
data_dir = '/home/rad/data/'+model+'/'+wind+'/'

sim =  caesar.load(data_dir+'Groups/'+model+'_'+snap+'.hdf5')

gal_sm = np.log10([i.masses['stellar'].in_units('Msun') for i in sim.central_galaxies])
gal_sfr = np.log10([i.sfr.in_units('Msun/yr') for i in sim.central_galaxies])

#mask = (gal_sm > sm_min) & (gal_sm < sm_max) & (gal_sfr > sfr_min) 
mask = []

for i in range(len(sim.central_galaxies)):
	mask.append(sfms_line(gal_sm[i], gal_sfr[i]))

