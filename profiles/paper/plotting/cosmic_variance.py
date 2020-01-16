import h5py
import numpy as np
import caesar
import sys
from plotting_methods import *

sys.path.append('/home/sapple/tools/')
from subdivide import octants, variance_jk

def cosmic_variance(profiles, pos, boxsize, quantity):
	octant_ids = octants(pos, boxsize)
	tukeys = np.zeros((8, len(profiles[0])))
	for i in range(8):
		i_using = np.concatenate(np.delete(octant_ids, i))
		tukeys[i], scale = tukey_biweight(profiles[i_using.astype('int')])
	mean_tukey = np.sum(tukeys, axis=0) / 8
   
	cosmic_var = variance_jk(tukeys, mean_tukey)
	cosmic_std = np.sqrt(cosmic_var)

	if quantity == 'sfr':
		mean_tukey[np.where(mean_tukey == 0.)[0]] = 1.e-6
	
	cosmic_std /= (np.log(10.)*mean_tukey)
	mean_tukey = np.log10(mean_tukey)
	
	return mean_tukey, cosmic_std

model = sys.argv[1]
wind = sys.argv[2]
snap = sys.argv[3]
sample_file = sys.argv[4]
angle = sys.argv[5]

with h5py.File(sample_file, 'r') as f:
    try:
        gals = f[model+'_'+snap][:]
    except KeyError:
        print('Need to identify galaxies first')

if 'gv' in sample_file.split('/', -1)[-1]:
		selection = 'gv'
elif 'sf' in sample_file.split('/', -1)[-1]:
		selection = 'sf'

basic_dir = '/home/sapple/simba_sizes/profiles/paper/with_dust_sizes/'
profs_file = basic_dir + model+'_'+snap+'/all_profiles_'+angle + '.h5'

results_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/data/'+selection
results_dir += selection
if angle == 'random_orientation':
	results_dir += '_rand'
elif angle == 'rotated_faceon':
	results_dir += '_rot'

mass_bins = [10., 10.5, 11.0]
bin_labels = ['10.0-10.5', '10.5-11.0', '>11.0']
mass_bins = [10., 10.6]
bin_labels = ['10.0-10.6','>10.6']

caesar_dir = '/home/rad/data/'+model+'/'+wind+'/Groups/'
sim =  caesar.load(caesar_dir+model+'_'+snap+'.hdf5', LoadHalo=False)
gal_cent = np.array([i.central for i in sim.galaxies])
gal_sm = np.log10(np.array([i.masses['stellar'].in_units('Msun') for i in sim.galaxies]))
gal_ids = np.array([False for i in range(len(sim.galaxies))])
gal_ids[gals] = True
gal_pos = np.array([i.pos.in_units('kpc/h') for i in sim.galaxies])
boxsize = sim.simulation.boxsize.in_units('kpc/h')

for i, b in enumerate(bin_labels):

	if i != len(mass_bins) -1:
			mass_mask = (gal_sm > mass_bins[i]) & (gal_sm < mass_bins[i+1])
	else:
			mass_mask = gal_sm > mass_bins[i]
	
	mask = mass_mask*gal_ids*gal_cent
	with h5py.File(profs_file, 'r') as f:
		cen_star_m = f['sm'].value[mask]
		cen_gas_sfr = f['gas_sfr'].value[mask]
		cen_gas_h1 = f['h1_m'].value[mask]
		cen_gas_h2 = f['h2_m'].value[mask]
		#cen_gal_ids = f['gal_ids'].value[mask]

	cen_gas_ssfr = cen_gas_sfr / cen_star_m
	cen_gas_fmol = cen_gas_h2 / cen_gas_h1
	cen_gas_sfe = cen_gas_sfr / cen_gas_h2
	cen_gas_fh2 = cen_gas_h2 / cen_star_m
	cen_pos = gal_pos[mask]

	mask = mass_mask*gal_ids*np.invert(gal_cent)
	with h5py.File(profs_file, 'r') as f:
		sat_star_m = f['sm'].value[mask]
		sat_gas_sfr = f['gas_sfr'].value[mask]
		sat_gas_h1 = f['h1_m'].value[mask]
		sat_gas_h2 = f['h2_m'].value[mask]
		#sat_gal_ids = f['gal_ids'].value[mask]

	sat_gas_ssfr = sat_gas_sfr / sat_star_m
	sat_gas_fmol = sat_gas_h2 / sat_gas_h1
	sat_gas_sfe = sat_gas_sfr / sat_gas_h2
	sat_gas_fh2 = sat_gas_h2 / sat_star_m
	sat_pos = gal_pos[mask]

	if i == 0:
		n = cen_star_m.shape[1]
		cen_ssfr_jk = np.zeros((len(bin_labels), n)); cen_ssfr_cv_err = np.zeros((len(bin_labels), n))
		cen_sfr_jk = np.zeros((len(bin_labels), n)); cen_sfr_cv_err = np.zeros((len(bin_labels), n))
		cen_h1_jk = np.zeros((len(bin_labels), n)); cen_h1_cv_err = np.zeros((len(bin_labels), n))
		cen_h2_jk = np.zeros((len(bin_labels), n)); cen_h2_cv_err = np.zeros((len(bin_labels), n))
		cen_fmol_jk = np.zeros((len(bin_labels), n)); cen_fmol_cv_err = np.zeros((len(bin_labels), n))
		cen_sfe_jk = np.zeros((len(bin_labels), n)); cen_sfe_cv_err = np.zeros((len(bin_labels), n))
		cen_fh2_jk = np.zeros((len(bin_labels), n)); cen_fh2_cv_err = np.zeros((len(bin_labels), n))


		sat_ssfr_jk = np.zeros((len(bin_labels), n)); sat_ssfr_cv_err = np.zeros((len(bin_labels), n))
		sat_sfr_jk = np.zeros((len(bin_labels), n)); sat_sfr_cv_err = np.zeros((len(bin_labels), n))
		sat_h1_jk = np.zeros((len(bin_labels), n)); sat_h1_cv_err = np.zeros((len(bin_labels), n))
		sat_h2_jk = np.zeros((len(bin_labels), n)); sat_h2_cv_err = np.zeros((len(bin_labels), n))
		sat_fmol_jk = np.zeros((len(bin_labels), n)); sat_fmol_cv_err = np.zeros((len(bin_labels), n))
		sat_sfe_jk = np.zeros((len(bin_labels), n)); sat_sfe_cv_err = np.zeros((len(bin_labels), n))
		sat_fh2_jk = np.zeros((len(bin_labels), n)); sat_fh2_cv_err = np.zeros((len(bin_labels), n))

		all_ssfr_jk = np.zeros((len(bin_labels), n)); all_ssfr_cv_err = np.zeros((len(bin_labels), n))
		all_sfr_jk = np.zeros((len(bin_labels), n)); all_sfr_cv_err = np.zeros((len(bin_labels), n))
		all_h1_jk = np.zeros((len(bin_labels), n)); all_h1_cv_err = np.zeros((len(bin_labels), n))
		all_h2_jk = np.zeros((len(bin_labels), n)); all_h2_cv_err = np.zeros((len(bin_labels), n))
		all_fmol_jk = np.zeros((len(bin_labels), n)); all_fmol_cv_err = np.zeros((len(bin_labels), n))
		all_sfe_jk = np.zeros((len(bin_labels), n)); all_sfe_cv_err = np.zeros((len(bin_labels), n))
		all_fh2_jk = np.zeros((len(bin_labels), n)); all_fh2_cv_err = np.zeros((len(bin_labels), n))

	cen_ssfr_jk[i], cen_ssfr_cv_err[i] = cosmic_variance(cen_gas_ssfr, cen_pos, boxsize, 'ssfr')  
	cen_sfr_jk[i], cen_sfr_cv_err[i] = cosmic_variance(cen_gas_sfr, cen_pos, boxsize, 'sfr')
	cen_h1_jk[i], cen_h1_cv_err[i] = cosmic_variance(cen_gas_h1, cen_pos, boxsize, 'h1')
	cen_h2_jk[i], cen_h2_cv_err[i] = cosmic_variance(cen_gas_h2, cen_pos, boxsize, 'h2')
	cen_fmol_jk[i], cen_fmol_cv_err[i] = cosmic_variance(cen_gas_fmol, cen_pos, boxsize, 'fmol')
	cen_sfe_jk[i], cen_sfe_cv_err[i] = cosmic_variance(cen_gas_sfe, cen_pos, boxsize, 'sfe')
	cen_fh2_jk[i], cen_fh2_cv_err[i] = cosmic_variance(cen_gas_fh2, cen_pos, boxsize, 'fh2')


	sat_ssfr_jk[i], sat_ssfr_cv_err[i] = cosmic_variance(sat_gas_ssfr, sat_pos, boxsize, 'ssfr')
	sat_sfr_jk[i], sat_sfr_cv_err[i] = cosmic_variance(sat_gas_sfr, sat_pos, boxsize, 'sfr')
	sat_h1_jk[i], sat_h1_cv_err[i] = cosmic_variance(sat_gas_h1, sat_pos, boxsize, 'h1')
	sat_h2_jk[i], sat_h2_cv_err[i] = cosmic_variance(sat_gas_h2, sat_pos, boxsize, 'h2')
	sat_fmol_jk[i], sat_fmol_cv_err[i] = cosmic_variance(sat_gas_fmol, sat_pos, boxsize, 'fmol')
	sat_sfe_jk[i], sat_sfe_cv_err[i] = cosmic_variance(sat_gas_sfe, sat_pos, boxsize, 'sfe')
	sat_fh2_jk[i], sat_fh2_cv_err[i] = cosmic_variance(sat_gas_fh2, sat_pos, boxsize, 'fh2')

	ssfr = np.concatenate((cen_gas_ssfr, sat_gas_ssfr))
	sfr = np.concatenate((cen_gas_sfr, sat_gas_sfr))
	h1 = np.concatenate((cen_gas_h1, sat_gas_h1))
	h2 = np.concatenate((cen_gas_h2, sat_gas_h2))
	fmol = np.concatenate((cen_gas_fmol, sat_gas_fmol))
	sfe = np.concatenate((cen_gas_sfe, sat_gas_sfe))
	fh2 = np.concatenate((cen_gas_fh2, sat_gas_fh2))
	all_pos = np.concatenate((cen_pos, sat_pos))

	all_ssfr_jk[i], all_ssfr_cv_err[i] = cosmic_variance(ssfr, all_pos, boxsize, 'ssfr')
	all_sfr_jk[i], all_sfr_cv_err[i] = cosmic_variance(sfr, all_pos, boxsize, 'sfr')
	all_h1_jk[i], all_h1_cv_err[i] = cosmic_variance(h1, all_pos, boxsize, 'h1')
	all_h2_jk[i], all_h2_cv_err[i] = cosmic_variance(h2, all_pos, boxsize, 'h2')
	all_fmol_jk[i], all_fmol_cv_err[i] = cosmic_variance(fmol, all_pos, boxsize, 'fmol')
	all_sfe_jk[i], all_sfe_cv_err[i] = cosmic_variance(sfe, all_pos, boxsize, 'sfe')
	all_fh2_jk[i], all_fh2_cv_err[i] = cosmic_variance(fh2, all_pos, boxsize, 'fh2')

	with h5py.File(results_dir+'_ssfr_data.h5', 'a') as f:
		f.create_dataset('cen_jk_'+bin_labels[i], data=np.array(cen_ssfr_jk[i]))
		f.create_dataset('cen_cv_err_'+bin_labels[i], data=np.array(cen_ssfr_cv_err[i]))
		f.create_dataset('sat_jk_'+bin_labels[i], data=np.array(sat_ssfr_jk[i]))
		f.create_dataset('sat_cv_err_'+bin_labels[i], data=np.array(sat_ssfr_cv_err[i]))
		f.create_dataset('all_jk_'+bin_labels[i], data=np.array(all_ssfr_jk[i]))
		f.create_dataset('all_cv_err_'+bin_labels[i], data=np.array(all_ssfr_cv_err[i]))

	with h5py.File(results_dir+'_sfr_data.h5', 'a') as f: 
		f.create_dataset('cen_jk_'+bin_labels[i], data=np.array(cen_sfr_jk[i]))
		f.create_dataset('cen_cv_err_'+bin_labels[i], data=np.array(cen_sfr_cv_err[i]))
		f.create_dataset('sat_jk_'+bin_labels[i], data=np.array(sat_sfr_jk[i]))
		f.create_dataset('sat_cv_err_'+bin_labels[i], data=np.array(sat_sfr_cv_err[i]))
		f.create_dataset('all_jk_'+bin_labels[i], data=np.array(all_sfr_jk[i]))
		f.create_dataset('all_cv_err_'+bin_labels[i], data=np.array(all_sfr_cv_err[i]))

	with h5py.File(results_dir+'_h1_data.h5', 'a') as f: 
		f.create_dataset('cen_jk_'+bin_labels[i], data=np.array(cen_h1_jk[i]))
		f.create_dataset('cen_cv_err_'+bin_labels[i], data=np.array(cen_h1_cv_err[i]))
		f.create_dataset('sat_jk_'+bin_labels[i], data=np.array(sat_h1_jk[i]))
		f.create_dataset('sat_cv_err_'+bin_labels[i], data=np.array(sat_h1_cv_err[i]))
		f.create_dataset('all_jk_'+bin_labels[i], data=np.array(all_h1_jk[i]))
		f.create_dataset('all_cv_err_'+bin_labels[i], data=np.array(all_h1_cv_err[i]))

	with h5py.File(results_dir+'_h2_data.h5', 'a') as f:
		f.create_dataset('cen_jk_'+bin_labels[i], data=np.array(cen_h2_jk[i]))
		f.create_dataset('cen_cv_err_'+bin_labels[i], data=np.array(cen_h2_cv_err[i]))
		f.create_dataset('sat_jk_'+bin_labels[i], data=np.array(sat_h2_jk[i]))
		f.create_dataset('sat_cv_err_'+bin_labels[i], data=np.array(sat_h2_cv_err[i]))
		f.create_dataset('all_jk_'+bin_labels[i], data=np.array(all_h2_jk[i]))
		f.create_dataset('all_cv_err_'+bin_labels[i], data=np.array(all_h2_cv_err[i]))

	with h5py.File(results_dir+'_fmol_data.h5', 'a') as f: 
		f.create_dataset('cen_jk_'+bin_labels[i], data=np.array(cen_fmol_jk[i]))
		f.create_dataset('cen_cv_err_'+bin_labels[i], data=np.array(cen_fmol_cv_err[i]))
		f.create_dataset('sat_jk_'+bin_labels[i], data=np.array(sat_fmol_jk[i]))
		f.create_dataset('sat_cv_err_'+bin_labels[i], data=np.array(sat_fmol_cv_err[i]))
		f.create_dataset('all_jk_'+bin_labels[i], data=np.array(all_fmol_jk[i]))
		f.create_dataset('all_cv_err_'+bin_labels[i], data=np.array(all_fmol_cv_err[i]))

	with h5py.File(results_dir+'_sfe_data.h5', 'a') as f: 
		f.create_dataset('cen_jk_'+bin_labels[i], data=np.array(cen_sfe_jk[i]))
		f.create_dataset('cen_cv_err_'+bin_labels[i], data=np.array(cen_sfe_cv_err[i]))
		f.create_dataset('sat_jk_'+bin_labels[i], data=np.array(sat_sfe_jk[i]))
		f.create_dataset('sat_cv_err_'+bin_labels[i], data=np.array(sat_sfe_cv_err[i]))
		f.create_dataset('all_jk_'+bin_labels[i], data=np.array(all_sfe_jk[i]))
		f.create_dataset('all_cv_err_'+bin_labels[i], data=np.array(all_sfe_cv_err[i]))

	with h5py.File(results_dir+'_fh2_data.h5', 'a') as f: 
		f.create_dataset('cen_jk_'+bin_labels[i], data=np.array(cen_fh2_jk[i]))
		f.create_dataset('cen_cv_err_'+bin_labels[i], data=np.array(cen_fh2_cv_err[i]))
		f.create_dataset('sat_jk_'+bin_labels[i], data=np.array(sat_fh2_jk[i]))
		f.create_dataset('sat_cv_err_'+bin_labels[i], data=np.array(sat_fh2_cv_err[i]))
		f.create_dataset('all_jk_'+bin_labels[i], data=np.array(all_fh2_jk[i]))
		f.create_dataset('all_cv_err_'+bin_labels[i], data=np.array(all_fh2_cv_err[i]))

"""
# for gals > 10.5:
bin_label = '>10.5'

with h5py.File(centrals_dir+'mask_3_all_profiles.h5', 'r') as f:
	cen_star_m = f['sm'].value
	cen_gas_sfr = f['gas_sfr'].value
	cen_gas_h1 = f['h1_m'].value
	cen_gas_h2 = f['h2_m'].value
	cen_gal_ids = f['gal_ids'].value

with h5py.File(centrals_dir+'mask_4_all_profiles.h5', 'r') as f:
	cen_star_m = np.concatenate((cen_star_m, f['sm'].value))
	cen_gas_sfr = np.concatenate((cen_gas_sfr, f['gas_sfr'].value))
	cen_gas_h1 = np.concatenate((cen_gas_h1,f['h1_m'].value))
	cen_gas_h2 = np.concatenate((cen_gas_h2,f['h2_m'].value))
	cen_gas_ids = np.concatenate((cen_gal_ids,f['gal_ids'].value))

with h5py.File(sats_dir+'mask_3_all_profiles.h5', 'r') as f:
	sat_star_m = f['sm'].value
	sat_gas_sfr = f['gas_sfr'].value
	sat_gas_h1 = f['h1_m'].value
	sat_gas_h2 = f['h2_m'].value
	sat_gal_ids = f['gal_ids'].value

with h5py.File(sats_dir+'mask_4_all_profiles.h5', 'r') as f:
	sat_star_m = np.concatenate((sat_star_m, f['sm'].value))
	sat_gas_sfr = np.concatenate((sat_gas_sfr, f['gas_sfr'].value))
	sat_gas_h1 = np.concatenate((sat_gas_h1, f['h1_m'].value))
	sat_gas_h2 = np.concatenate((sat_gas_h2, f['h2_m'].value))
	sat_gas_ids = np.concatenate((sat_gal_ids,f['gal_ids'].value))

# centrals:
cen_no_gals = len(cen_gas_sfr)
cen_gas_sfe = cen_gas_sfr / cen_gas_h2
cen_gas_fh2 = cen_gas_h2 / cen_star_m
cen_gas_ssfr = cen_gas_sfr / cen_star_m
cen_gas_fmol = cen_gas_h2 / cen_gas_h1
cen_pos = gal_pos[cen_gal_ids]

cen_ssfr_jk, cen_ssfr_cv_err = cosmic_variance(cen_gas_ssfr, cen_pos, boxsize, 'ssfr')
cen_sfr_jk, cen_sfr_cv_err = cosmic_variance(cen_gas_sfr, cen_pos, boxsize, 'sfr')
cen_h1_jk, cen_h1_cv_err = cosmic_variance(cen_gas_h1, cen_pos, boxsize, 'h1')
cen_h2_jk, cen_h2_cv_err = cosmic_variance(cen_gas_h2, cen_pos, boxsize, 'h2')
cen_fmol_jk, cen_fmol_cv_err = cosmic_variance(cen_gas_fmol, cen_pos, boxsize, 'fmol')
cen_sfe_jk, cen_sfe_cv_err = cosmic_variance(cen_gas_sfe, cen_pos, boxsize, 'sfe')
cen_fh2_jk, cen_fh2_cv_err = cosmic_variance(cen_gas_fh2, cen_pos, boxsize, 'fh2')

# satellites:
sat_no_gals = len(sat_gas_sfr)
sat_gas_sfe = sat_gas_sfr / sat_gas_h2
sat_gas_fh2 = sat_gas_h2 / sat_star_m
sat_gas_ssfr = sat_gas_sfr / sat_star_m
sat_gas_fmol = sat_gas_h2 / sat_gas_h1
sat_pos = gal_pos[sat_gal_ids]

sat_ssfr_jk, sat_ssfr_cv_err = cosmic_variance(sat_gas_ssfr, sat_pos, boxsize, 'ssfr')
sat_sfr_jk, sat_sfr_cv_err = cosmic_variance(sat_gas_sfr, sat_pos, boxsize, 'sfr')
sat_h1_jk, sat_h1_cv_err = cosmic_variance(sat_gas_h1, sat_pos, boxsize, 'h1')
sat_h2_jk, sat_h2_cv_err = cosmic_variance(sat_gas_h2, sat_pos, boxsize, 'h2')
sat_fmol_jk, sat_fmol_cv_err = cosmic_variance(sat_gas_fmol, sat_pos, boxsize, 'fmol')
sat_sfe_jk, sat_sfe_cv_err = cosmic_variance(sat_gas_sfe, sat_pos, boxsize, 'sfe')
sat_fh2_jk, sat_fh2_cv_err = cosmic_variance(sat_gas_fh2, sat_pos, boxsize, 'fh2')

with h5py.File(results_dir+'_sfr_data.h5', 'a') as f:
		#f.create_dataset('cen_no_gals_'+bin_label, data=np.array(cen_no_gals))
		#f.create_dataset('sat_no_gals_'+bin_label, data=np.array(sat_no_gals))
		f.create_dataset('cen_jk_'+bin_label, data=np.array(cen_sfr_jk))
		f.create_dataset('cen_cv_err_'+bin_label, data=np.array(cen_sfr_cv_err))
		f.create_dataset('sat_jk_'+bin_label, data=np.array(sat_sfr_jk))
		f.create_dataset('sat_cv_err_'+bin_label, data=np.array(sat_sfr_cv_err))

with h5py.File(results_dir+'_ssfr_data.h5', 'a') as f:
		#f.create_dataset('cen_no_gals_'+bin_label, data=np.array(cen_no_gals))
		#f.create_dataset('sat_no_gals_'+bin_label, data=np.array(sat_no_gals))
		f.create_dataset('cen_jk_'+bin_label, data=np.array(cen_ssfr_jk))
		f.create_dataset('cen_cv_err_'+bin_label, data=np.array(cen_ssfr_cv_err))
		f.create_dataset('sat_jk_'+bin_label, data=np.array(sat_ssfr_jk))
		f.create_dataset('sat_cv_err_'+bin_label, data=np.array(sat_ssfr_cv_err))

with h5py.File(results_dir+'_h1_data.h5', 'a') as f:
		#f.create_dataset('cen_no_gals_'+bin_label, data=np.array(cen_no_gals))
		#f.create_dataset('sat_no_gals_'+bin_label, data=np.array(sat_no_gals))
		f.create_dataset('cen_jk_'+bin_label, data=np.array(cen_h1_jk))
		f.create_dataset('cen_cv_err_'+bin_label, data=np.array(cen_h1_cv_err))
		f.create_dataset('sat_jk_'+bin_label, data=np.array(sat_h1_jk))
		f.create_dataset('sat_cv_err_'+bin_label, data=np.array(sat_h1_cv_err))

with h5py.File(results_dir+'_h2_data.h5', 'a') as f:
		#f.create_dataset('cen_no_gals_'+bin_label, data=np.array(cen_no_gals))
		#f.create_dataset('sat_no_gals_'+bin_label, data=np.array(sat_no_gals))
		f.create_dataset('cen_jk_'+bin_label, data=np.array(cen_h2_jk))
		f.create_dataset('cen_cv_err_'+bin_label, data=np.array(cen_h2_cv_err))
		f.create_dataset('sat_jk_'+bin_label, data=np.array(sat_h2_jk))
		f.create_dataset('sat_cv_err_'+bin_label, data=np.array(sat_h2_cv_err))

with h5py.File(results_dir+'_fmol_data.h5', 'a') as f:
		#f.create_dataset('cen_no_gals_'+bin_label, data=np.array(cen_no_gals))
		#f.create_dataset('sat_no_gals_'+bin_label, data=np.array(sat_no_gals))
		f.create_dataset('cen_jk_'+bin_label, data=np.array(cen_fmol_jk))
		f.create_dataset('cen_cv_err_'+bin_label, data=np.array(cen_fmol_cv_err))
		f.create_dataset('sat_jk_'+bin_label, data=np.array(sat_fmol_jk))
		f.create_dataset('sat_cv_err_'+bin_label, data=np.array(sat_fmol_cv_err))

with h5py.File(results_dir+'_sfe_data.h5', 'a') as f:
		#f.create_dataset('cen_no_gals_'+bin_label, data=np.array(cen_no_gals))
		#f.create_dataset('sat_no_gals_'+bin_label, data=np.array(sat_no_gals))
		f.create_dataset('cen_jk_'+bin_label, data=np.array(cen_sfe_jk))
		f.create_dataset('cen_cv_err_'+bin_label, data=np.array(cen_sfe_cv_err))
		f.create_dataset('sat_jk_'+bin_label, data=np.array(sat_sfe_jk))
		f.create_dataset('sat_cv_err_'+bin_label, data=np.array(sat_sfe_cv_err))

with h5py.File(results_dir+'_fh2_data.h5', 'a') as f:
		#f.create_dataset('cen_no_gals_'+bin_label, data=np.array(cen_no_gals))
		#f.create_dataset('sat_no_gals_'+bin_label, data=np.array(sat_no_gals))
		f.create_dataset('cen_jk_'+bin_label, data=np.array(cen_fh2_jk))
		f.create_dataset('cen_cv_err_'+bin_label, data=np.array(cen_fh2_cv_err))
		f.create_dataset('sat_jk_'+bin_label, data=np.array(sat_fh2_jk))
		f.create_dataset('sat_cv_err_'+bin_label, data=np.array(sat_fh2_cv_err))

"""
