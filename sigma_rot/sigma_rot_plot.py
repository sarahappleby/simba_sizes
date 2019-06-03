import matplotlib.pyplot as plt
import h5py
import numpy as np

if __name__ == '__main__':

	model = 'm50n512'
	wind = 's50j7k'
	snap = '151'

	results_dir = '/home/sapple/simba_sizes/sigma_rot/'

	with h5py.File(results_dir+model+'_'+snap+'_kinematics.h5', 'r') as f:
		sigv_h1 = f['sigmav_h1'].value
		sigv_h2 = f['sigmav_h2'].value
		#vrot_gravity_h1 = f['vrot_gravity_h1'].value # these give nans
		#vrot_gravity_h2 = f['vrot_gravity_h2'].value # these give nans
		vrot_h1 = f['vrot_h1'].value
		vrot_h2 = f['vrot_h2'].value
		gal_ids = f['gal_ids'].value

	with h5py.File(results_dir+model+'_'+snap+'_caesar_quantities.h5', 'r') as f:
		gal_sm = f['stellar_mass'].value[gal_ids]
		gal_bm = f['bh_mass'].value[gal_ids]
		gal_sfr = f['sfr'].value[gal_ids]
	
	vsigma_h1 = vrot_h1 / sigv_h1
	vsigma_h2 = vrot_h2 / sigv_h2
	ssfr = gal_sfr / gal_sm
	c = np.log10(ssfr /1.e-9)

	plt.scatter(np.log10(gal_sm), vsigma_h1, c=c, cmap='jet_r', s=2)
	plt.colorbar(label='log sSFR (Gyr)')
	plt.clim(-2, )
	plt.xlabel('M*')
	plt.ylabel('v/sigma')
	plt.show()