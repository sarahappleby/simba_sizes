import h5py
import numpy as np
import sys
from plotting_methods import *

model = sys.argv[1]
wind = sys.argv[2]
snap = sys.argv[3]
selection = sys.argv[4]
gals = sys.argv[5]
angle = sys.argv[6]

if selection =='gv':
    name = 'green_valley'
elif selection == 'sf':
    name = 'star_forming'
basic_dir = '/home/sapple/simba_sizes/profiles/paper/'

centrals_dir = basic_dir + 'centrals/'+model+'_'+snap+'/'+wind+'/'+name+'/'+angle+'/'
sats_dir = basic_dir + 'satellites/'+model+'_'+snap+'/'+wind+'/'+name+'/'+angle+'/'
results_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/data/'

results_dir += selection

if 'random' in centrals_dir:
    results_dir += '_rand'
elif 'rotated' in centrals_dir:
    results_dir += '_rot'

masks = [2, 3, 4]
bin_labels = ['10.0-10.5', '10.5-11.0', '>11.0']

cen_no_gals = np.zeros(len(masks))
sat_no_gals = np.zeros(len(masks))
all_no_gals = np.zeros(len(masks))
for i, m in enumerate(masks):

    with h5py.File(centrals_dir+'mask_'+str(m)+'_all_profiles.h5', 'r') as f:
        cen_gas_h2_z = f['h2-weighted_z'].value

    if gals == 'all':
        with h5py.File(sats_dir+'mask_'+str(m)+'_all_profiles.h5', 'r') as f:
            sat_gas_h2_z = f['h2-weighted_z'].value
        

    if i == 0:
        n = cen_gas_h2_z.shape[1]
        cen_h2_z_tukey = np.zeros((len(masks), n)); cen_h2_z_large_scale = np.zeros((len(masks), n)); cen_h2_z_small_scale = np.zeros((len(masks), n))

        sat_h2_z_tukey = np.zeros((len(masks), n)); sat_h2_z_large_scale = np.zeros((len(masks), n)); sat_h2_z_small_scale = np.zeros((len(masks), n))

        if gals == 'all':
            all_h2_z_tukey = np.zeros((len(masks), n)); all_h2_z_large_scale = np.zeros((len(masks), n)); all_h2_z_small_scale = np.zeros((len(masks), n))

    # centrals:
    cen_no_gals[i] = len(cen_gas_h2_z)

    tukey, scale = tukey_biweight(cen_gas_h2_z)
    cen_h2_z_tukey[i] = np.log10(tukey)
    cen_h2_z_large_scale[i] = scale / (np.log(10.)*tukey)
    cen_h2_z_small_scale[i] = scale / (np.sqrt(cen_no_gals[i])* np.log(10.)*tukey)

    if gals == 'all':
        # satellites:
        sat_no_gals[i] = len(sat_gas_h2_z)

        tukey, scale = tukey_biweight(sat_gas_h2_z)
        sat_h2_z_tukey[i] = np.log10(tukey)
        sat_h2_z_large_scale[i] = scale / (np.log(10.)*tukey)
        sat_h2_z_small_scale[i] = scale / (np.sqrt(cen_no_gals[i])* np.log(10.)*tukey)

        h2_z = np.concatenate((cen_gas_h2_z, sat_gas_h2_z))
        all_no_gals[i] = len(h2_z)

        tukey, scale = tukey_biweight(h2_z)
        all_h2_z_tukey[i] = np.log10(tukey)
        all_h2_z_large_scale[i] = scale / (np.log(10.)*tukey)
        all_h2_z_small_scale[i] = scale / (np.sqrt(all_no_gals[i])* np.log(10.)*tukey)


    with h5py.File(results_dir+'_h2-weighted_z_data.h5', 'a') as f:
        f.create_dataset('cen_no_gals_'+bin_labels[i], data=np.array(cen_no_gals[i]))
        f.create_dataset('cen_tukey_'+bin_labels[i], data=np.array(cen_h2_z_tukey[i]))
        f.create_dataset('cen_large_scale_'+bin_labels[i], data=np.array(cen_h2_z_large_scale[i]))
        f.create_dataset('cen_small_scale_'+bin_labels[i], data=np.array(cen_h2_z_small_scale[i]))
        if gals == 'all':
            f.create_dataset('all_no_gals_'+bin_labels[i], data=np.array(all_no_gals[i]))
            f.create_dataset('all_tukey_'+bin_labels[i], data=np.array(all_h2_z_tukey[i]))
            f.create_dataset('all_large_scale_'+bin_labels[i], data=np.array(all_h2_z_large_scale[i]))
            f.create_dataset('all_small_scale_'+bin_labels[i], data=np.array(all_h2_z_small_scale[i]))
            f.create_dataset('sat_no_gals_'+bin_labels[i], data=np.array(sat_no_gals[i]))
            f.create_dataset('sat_tukey_'+bin_labels[i], data=np.array(sat_h2_z_tukey[i]))
            f.create_dataset('sat_large_scale_'+bin_labels[i], data=np.array(sat_h2_z_large_scale[i]))
            f.create_dataset('sat_small_scale_'+bin_labels[i], data=np.array(sat_h2_z_small_scale[i]))
"""
# for gals > 10.5:
bin_label = '>10.5'

with h5py.File(centrals_dir+'mask_3_all_profiles.h5', 'r') as f:
    cen_gas_h2_z = f['h2_z'].value

with h5py.File(centrals_dir+'mask_4_all_profiles.h5', 'r') as f:
    cen_gas_h2_z = np.concatenate((cen_gas_h2_z,f['h2_z'].value))

with h5py.File(sats_dir+'mask_3_all_profiles.h5', 'r') as f:
    sat_gas_h2_z = f['h2_z'].value

with h5py.File(sats_dir+'mask_4_all_profiles.h5', 'r') as f:
    sat_gas_h2_z = np.concatenate((sat_gas_h2_z, f['h2_z'].value))

# centrals:
cen_no_gals = len(cen_gas_h2_z)

tukey, scale = tukey_biweight(cen_gas_h2_z)
cen_h2_z_tukey = np.log10(tukey)
cen_h2_z_large_scale = scale / (np.log(10.)*tukey)
cen_h2_z_small_scale = scale / (np.sqrt(cen_no_gals)* np.log(10.)*tukey)

# satellites:
sat_no_gals = len(sat_gas_h2_z)

tukey, scale = tukey_biweight(sat_gas_h2_z)
sat_h2_z_tukey = np.log10(tukey)
sat_h2_z_large_scale = scale / (np.log(10.)*tukey)
sat_h2_z_small_scale = scale / (np.sqrt(sat_no_gals)* np.log(10.)*tukey)

with h5py.File(results_dir+'_h2_z_data.h5', 'a') as f:
        f.create_dataset('cen_no_gals_'+bin_label, data=np.array(cen_no_gals))
        f.create_dataset('sat_no_gals_'+bin_label, data=np.array(sat_no_gals))
        f.create_dataset('cen_tukey_'+bin_label, data=np.array(cen_h2_z_tukey))
        f.create_dataset('sat_tukey_'+bin_label, data=np.array(sat_h2_z_tukey))
        f.create_dataset('cen_large_scale_'+bin_label, data=np.array(cen_h2_z_large_scale))
        f.create_dataset('sat_large_scale_'+bin_label, data=np.array(sat_h2_z_large_scale))
        f.create_dataset('cen_small_scale_'+bin_label, data=np.array(cen_h2_z_small_scale))
        f.create_dataset('sat_small_scale_'+bin_label, data=np.array(sat_h2_z_small_scale))
"""
