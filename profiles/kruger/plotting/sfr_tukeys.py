import h5py
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
from plotting_methods import tukey_biweight

profile_dir = sys.argv[1]
results_dir = profile_dir

if 'satellites' in profile_dir:
    satellites = True
    masks = [2, 3]
    bin_labels = ['10.0-10.5', '>10.5']
else:
    satellites = False
    masks = [2, 3, 4]
    bin_labels = ['10.0-10.5', '10.5-11.0', '>11.0']

no_gals = np.zeros(len(bin_labels))

for i, m in enumerate(masks):
    with h5py.File(profile_dir+'mask_'+str(m)+'_all_profiles.h5', 'r') as f:
        sfr = f['gas_sfr'].value
    if satellites & (m == 3):
        with h5py.File(profile_dir+'mask_'+str(m+1)+'_all_profiles.h5', 'r') as f:
            ages = np.concatenate((sfr, f['gas_sfr'].value))

    if i == 0:
        n = sfr.shape[1]
        sfr_tukey = np.zeros((len(bin_labels), n)); sfr_large_scale = np.zeros((len(bin_labels), n)); sfr_small_scale = np.zeros((len(bin_labels), n))
    no_gals[i] = len(sfr)

    tukey, scale = tukey_biweight(sfr)
    sfr_tukey[i] = np.log10(tukey)
    sfr_large_scale[i] = scale / (np.log(10.)*tukey)
    sfr_small_scale[i] = scale / (np.sqrt(no_gals[i])* np.log(10.)*tukey)

    with h5py.File(results_dir+'sfr_data.h5', 'a') as f:
        f.create_dataset('no_gals_'+bin_labels[i], data=np.array(no_gals[i]))
        f.create_dataset('tukey_'+bin_labels[i], data=np.array(sfr_tukey[i]))
        f.create_dataset('large_scale_'+bin_labels[i], data=np.array(sfr_large_scale[i]))
        f.create_dataset('small_scale_'+bin_labels[i], data=np.array(sfr_small_scale[i]))

