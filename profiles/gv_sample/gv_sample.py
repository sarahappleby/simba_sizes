import caesar
import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
import os

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

gv_ssfr_min = -13.

masses = [10., 10.5, 11.]


model = sys.argv[1]
wind = sys.argv[2]
snap = sys.argv[3]

h5_dir = '/home/sapple/simba_sizes/profiles/gv_sample/'+wind + '/'
if not os.path.exists(h5_dir):
	os.makedirs(h5_dir)
results_dir = h5_dir +model+'_'+snap+'/'
if not os.path.exists(results_dir):
	os.makedirs(results_dir)

data_dir = '/home/rad/data/'+model+'/'+wind+'/'

if model == 'm50n512' and wind == 's50j7k':
	sim =  caesar.load(data_dir+'Groups/caesar_old/'+model+'_'+snap+'.hdf5', LoadHalo=False)
else:
	sim =  caesar.load(data_dir+'Groups/'+model+'_'+snap+'.hdf5', LoadHalo=False)

gal_cent = np.array([i.central for i in sim.galaxies])
gal_sm = np.log10([i.masses['stellar'].in_units('Msun') for i in sim.galaxies])
gal_sfr = np.log10([i.sfr.in_units('Msun/yr') for i in sim.galaxies])
gal_ssfr = gal_sfr / gal_sm

"""
Get sample of the green valley from Belfiore 18
"""

ssfr_mask = (gal_ssfr > gv_ssfr_min) 

sfms_mask = []	
for i in range(len(sim.galaxies)):
	sfms_mask.append(sfms_line_check(gal_sm[i], gal_sfr[i]))

sm_plot = np.arange(9.0, 12., 0.5)
sfr_plot = sfms_line(sm_plot)

gv_mask = sfms_mask*ssfr_mask

plt.plot(sm_plot, sfr_plot, ls='--', lw=1.5, c='k', label='Belfiore 2018')
plt.scatter(gal_sm[gv_mask], gal_sfr[gv_mask], s=1, c='g', label='Green Valley')
plt.scatter(gal_sm[np.invert(sfms_mask)], gal_sfr[np.invert(sfms_mask)], s=1, c='b', label='SFMS')
plt.xlim(9., 11.5)
plt.legend()
plt.xlabel('log M*')
plt.ylabel('log SFR')
plt.savefig(results_dir+'b18_sample_all.png')
plt.clf()

cent_mask = gv_mask*gal_cent
plt.plot(sm_plot, sfr_plot, ls='--', lw=1.5, c='k', label='Belfiore 2018')
plt.scatter(gal_sm[cent_mask], gal_sfr[cent_mask], s=1, c='g', label='Green Valley')
plt.scatter(gal_sm[np.invert(sfms_mask)*gal_cent], gal_sfr[np.invert(sfms_mask)*gal_cent], s=1, c='b', label='SFMS')
plt.xlim(9., 11.5)
plt.legend()
plt.xlabel('log M*')
plt.ylabel('log SFR')
plt.savefig(results_dir+'b18_sample_centrals.png')
plt.clf()

for i in range(len(masses)):
	if i != 2:
		sm_mask = (gal_sm > masses[i]) & (gal_sm < masses[i+1])
	else:
		sm_mask = gal_sm > masses[i]
	mask = sfms_mask*ssfr_mask*sm_mask

	plt.plot(sm_plot, sfr_plot, ls='--', lw=1.5, c='k', label='Belfiore 2018')
	plt.scatter(gal_sm[mask], gal_sfr[mask], s=1, c='g', label='Green Valley')
	plt.scatter(gal_sm[np.invert(sfms_mask)*sm_mask], gal_sfr[np.invert(sfms_mask)*sm_mask], s=1, c='b', label='SFMS')
	plt.xlim(9., 11.5)
	plt.legend()
	plt.xlabel('log M*')
	plt.ylabel('log SFR')
	plt.savefig(results_dir+'b18_sample_mask_'+str(i)+'.png')
	plt.clf()

	mask *= gal_cent
	plt.plot(sm_plot, sfr_plot, ls='--', lw=1.5, c='k', label='Belfiore 2018')
	plt.scatter(gal_sm[mask], gal_sfr[mask], s=1, c='g', label='Green Valley')
	plt.scatter(gal_sm[np.invert(sfms_mask)*sm_mask*gal_cent], gal_sfr[np.invert(sfms_mask)*sm_mask*gal_cent], s=1, c='b', label='SFMS')
	plt.xlim(9., 11.5)
	plt.legend()
	plt.xlabel('log M*')
	plt.ylabel('log SFR')
	plt.savefig(results_dir+'b18_sample_centrals_mask_'+str(i)+'.png')
	plt.clf()

with h5py.File(h5_dir+'gv_samples.h5', 'a') as f:
	f.create_dataset(model+'_'+snap, data=np.array(gv_mask))
