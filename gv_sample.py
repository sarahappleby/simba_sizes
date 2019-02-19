import caesar
import numpy as np
import matplotlib.pyplot as plt
import h5py

#from matplotlib import rc
#rc('text', usetex=True)

def sfms_line(x0, a=0.73, b=-7.7):
	return x0*a + b

def sfms_line_check(x0, y0):
	y = sfms_line(x0)
	if y > y0:
		return True
	else:
		return False

# in belfiore 18, fig 2:
sfr_min = -1.5
sm_min = 10.

model = 'm50n512'
wind = 's50j7k'
snap = '151'

results_dir = '/home/sapple/simba_sizes/'
data_dir = '/home/rad/data/'+model+'/'+wind+'/'

sim =  caesar.load(data_dir+'Groups/'+model+'_'+snap+'.hdf5')

gal_sm = np.log10([i.masses['stellar'].in_units('Msun') for i in sim.central_galaxies])
gal_sfr = np.log10([i.sfr.in_units('Msun/yr') for i in sim.central_galaxies])

"""
Get sample of the green valley from Belfiore 18
"""

len_mask = np.array([len(i.glist) > 256 for i in sim.central_galaxies])
sfr_mask = (gal_sfr > sfr_min) 
sfms_mask = []
	
for i in range(len(sim.central_galaxies)):
	sfms_mask.append(sfms_line_check(gal_sm[i], gal_sfr[i]))

sm_plot = np.arange(9.0, 12., 0.5)
sfr_plot = sfms_line(sm_plot)

mask = sfms_mask*sfr_mask
plt.plot(sm_plot, sfr_plot, ls='--', lw=1.5, c='k', label='Belfiore 2018')
plt.scatter(gal_sm[mask], gal_sfr[mask], s=1, c='g', label='Green Valley')
plt.scatter(gal_sm[np.invert(sfms_mask)], gal_sfr[np.invert(sfms_mask)], s=1, c='b', label='SFMS')
plt.xlim(9., 11.5)
plt.legend()
plt.xlabel('log M*')
plt.ylabel('log SFR')
plt.savefig(results_dir+'b18_sample.png')
plt.clf()

mask *= len_mask
plt.plot(sm_plot, sfr_plot, ls='--', lw=1.5, c='k', label='Belfiore 2018')
plt.scatter(gal_sm[mask], gal_sfr[mask], s=1, c='g', label='Green Valley')
plt.scatter(gal_sm[np.invert(sfms_mask)], gal_sfr[np.invert(sfms_mask)], s=1, c='b', label='SFMS')
plt.xlim(9., 11.5)
plt.legend()
plt.xlabel('log M*')
plt.ylabel('log SFR')
plt.savefig(results_dir+'b18_sample_glist.png')
plt.clf()

"""
Get sample just from region of M*/SFR space and len mask
"""

ms_mask = (gal_sm > sm_min)
mask = len_mask*ms_mask*sfr_mask
plt.plot(sm_plot, sfr_plot, ls='--', lw=1.5, c='k', label='Belfiore 2018')
plt.scatter(gal_sm[mask], gal_sfr[mask], s=1, c='g', label='Green Valley')
plt.xlim(9., )
plt.legend()
plt.xlabel('log M*')
plt.ylabel('log SFR')
plt.savefig(results_dir+'basic_sample.png')
plt.clf()


with h5py.File(results_dir+'samples.h5', 'a') as f:
	f.create_dataset('len_mask', data=np.array(len_mask))
	f.create_dataset('ms_mask', data=np.array(sm_mask))
	f.create_dataset('sfr_mask', data=np.array(sfr_mask))
	f.create_dataset('sfms_mask', data=np.array(sfms_mask))