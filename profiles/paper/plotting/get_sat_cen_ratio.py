import h5py
import numpy as np
import sys
import caesar
from plotting_methods import *
sys.path.append('/home/sapple/tools/')
from subdivide import octants, variance_jk


def cosmic_variance_nolog(profiles, pos, boxsize, quantity):
	octant_ids = octants(pos, boxsize)
	tukeys = np.zeros((8, len(profiles[0])))
	for i in range(8):
		i_using = np.concatenate(np.delete(octant_ids, i))
		tukeys[i], scale = tukey_biweight(profiles[i_using.astype('int')])
	mean_tukey = np.sum(tukeys, axis=0) / 8

	cosmic_var = variance_jk(tukeys, mean_tukey)
	cosmic_std = np.sqrt(cosmic_var)

	return mean_tukey, cosmic_std


model = 'm100n1024'
wind = 's50j7k'
snap = '151'
sample_file = '../../gv_sample/belfiore/s50j7k/selection_2/sf_samples.h5'

with h5py.File(sample_file, 'r') as f:
	try:
		gals = f[model+'_'+snap].value
	except KeyError:
		print('Need to identify galaxies first')

if 'gv' in sample_file.split('/', -1)[-1]:
	selection = 'gv'
	name = 'green_valley'
elif 'sf' in sample_file.split('/', -1)[-1]:
	selection = 'sf'
	name = 'star_forming'

mass_bins = [10., 10.6,]
bin_labels = ['10.0-10.6', '>10.6']
angle = 'rotated_faceon'

masks = bin_labels

basic_dir = '/home/sapple/simba_sizes/profiles/paper/high_redshift/halfradius_units/'
basic_dir = '/home/sapple/simba_sizes/profiles/paper/high_redshift/physical_units/'
basic_dir = '/home/sapple/simba_sizes/profiles/paper/'

profs_file = basic_dir + 'all_profs/'+name+'_'+angle + '.h5'
results_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/data/'

results_dir += selection + '_rot' 

data_dir = '/home/rad/data/'+model+'/'+wind+'/'
sim =  caesar.load(data_dir+'Groups/'+model+'_'+snap+'.hdf5', LoadHalo=False)
gal_sm = np.log10(np.array([i.masses['stellar'].in_units('Msun') for i in sim.galaxies]))
gal_cent = np.array([i.central for i in sim.galaxies])
gal_ids = np.array([False for i in range(len(sim.galaxies))])
gal_ids[gals] = True

gal_pos = np.array([i.pos.in_units('kpc/h') for i in sim.galaxies])
boxsize = sim.simulation.boxsize.in_units('kpc/h')

cen_no_gals = np.zeros(len(bin_labels))
sat_no_gals = np.zeros(len(bin_labels))

for i, b in enumerate(bin_labels):

	if i != len(mass_bins) -1:
		mass_mask = (gal_sm > mass_bins[i]) & (gal_sm < mass_bins[i+1])
	else:
		mass_mask = gal_sm > mass_bins[i]

	mask = mass_mask*gal_ids*gal_cent
	with h5py.File(profs_file, 'r') as f:
		cen_star_m = f['sm'].value[mask]
		cen_gas_sfr = f['sfr'].value[mask]
		cen_gas_h1 = f['h1_m'].value[mask]
	cen_pos = gal_pos[mask]

	mask = mass_mask*gal_ids*np.invert(gal_cent)
	with h5py.File(profs_file, 'r') as f:
		sat_star_m = f['sm'].value[mask]
		sat_gas_sfr = f['sfr'].value[mask]
		sat_gas_h1 = f['h1_m'].value[mask]
	sat_pos = gal_pos[mask]
	
	cen_no_gals[i] = len(cen_gas_sfr)
	cen_gas_ssfr = cen_gas_sfr / cen_star_m
	sat_no_gals[i] = len(sat_gas_sfr)
	sat_gas_ssfr = sat_gas_sfr / sat_star_m

	#ssfr
	cen_tukey, cen_scale = tukey_biweight(cen_gas_ssfr)
	cen_ssfr_small_scale = cen_scale / np.sqrt(cen_no_gals[i])
	cen_ssfr_jk, cen_ssfr_cv_err = cosmic_variance_nolog(cen_gas_ssfr, cen_pos, boxsize, 'ssfr')
	cen_ssfr_err = np.sqrt(cen_ssfr_cv_err**2 + cen_ssfr_small_scale**2)


	sat_tukey, sat_scale = tukey_biweight(sat_gas_ssfr)
	sat_ssfr_small_scale = sat_scale / np.sqrt(sat_no_gals[i])
	sat_ssfr_jk, sat_ssfr_cv_err = cosmic_variance_nolog(sat_gas_ssfr, sat_pos, boxsize, 'ssfr')
	sat_ssfr_err = np.sqrt(sat_ssfr_cv_err**2 + sat_ssfr_small_scale**2)
	
	ssfr_ratio = np.log10(sat_tukey) - np.log10(cen_tukey)
	sat_log_ssfr_err = sat_ssfr_err / (sat_tukey * np.log(10.))
	cen_log_ssfr_err = cen_ssfr_err / (cen_tukey * np.log(10.))
	ssfr_ratio_err = np.sqrt(sat_log_ssfr_err**2 + cen_log_ssfr_err**2)

	#h1
	cen_tukey, cen_scale = tukey_biweight(cen_gas_h1)
	cen_h1_small_scale = cen_scale / np.sqrt(cen_no_gals[i])
	cen_h1_jk, cen_h1_cv_err = cosmic_variance_nolog(cen_gas_h1, cen_pos, boxsize, 'h1')
	cen_h1_err = np.sqrt(cen_h1_cv_err**2 + cen_h1_small_scale**2)

	sat_tukey, sat_scale = tukey_biweight(sat_gas_h1)
	sat_h1_small_scale = sat_scale / np.sqrt(sat_no_gals[i])
	sat_h1_jk, sat_h1_cv_err = cosmic_variance_nolog(sat_gas_h1, sat_pos, boxsize, 'h1')
	sat_h1_err = np.sqrt(sat_h1_cv_err**2 + sat_h1_small_scale**2)
		
	h1_ratio = np.log10(sat_tukey) - np.log10(cen_tukey)
	sat_log_h1_err = sat_h1_err / (sat_tukey * np.log(10.))
	cen_log_h1_err = cen_h1_err / (cen_tukey * np.log(10.))
	h1_ratio_err = np.sqrt(sat_log_h1_err**2 + cen_log_h1_err**2)



	with h5py.File(results_dir+'_sat_cent_ratio_data.h5', 'a') as f:
			f.create_dataset('cen_no_gals_'+b, data=np.array(cen_no_gals[i]))
			f.create_dataset('sat_no_gals_'+b, data=np.array(sat_no_gals[i]))
			f.create_dataset('ssfr_ratio_'+b, data=np.array(ssfr_ratio))
			f.create_dataset('ssfr_err_'+b, data=np.array(ssfr_ratio_err))
			f.create_dataset('h1_ratio_'+b, data=np.array(h1_ratio))
			f.create_dataset('h1_err_'+b, data=np.array(h1_ratio_err))

