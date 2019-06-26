import h5py
import numpy as np 
import matplotlib.pyplot as plt
import sys
import os
from plotting_methods import tukey_biweight

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 12})

profile_dir = sys.argv[1]
results_dir = profile_dir

if 'satellites' in profile_dir: 
    satellites = True
    masks = [2, 3]
    bin_labels = ['10.0-10.5', '>10.5']
    colors = ['b', 'm']
else: 
    satellites = False
    masks = [2, 3, 4]
    bin_labels = ['10.0-10.5', '10.5-11.0', '>11.0']
    colors = ['b', 'm', 'r']

no_gals = np.zeros(len(bin_labels))

for i, m in enumerate(masks):

	with h5py.File(profile_dir+'mask_'+str(m)+'_all_profiles.h5', 'r') as f:
    		star_m = f['sm'].value
                gas_m = f['gm'].value
		gas_h1 = f['h1'].value
		gas_h2 = f['h2'].value
        
        if satellites & (m == 3):
                with h5py.File(profile_dir+'mask_'+str(m+1)+'_all_profiles.h5', 'r') as f:
                    star_m = np.concatenate((star_m, f['sm'].value))
                    gas_m = np.concatenate((gas_m, f['gm'].value))
                    gas_h1 = np.concatenate((gas_h1, f['h1'].value))
                    gas_h2 = np.concatenate((gas_h2, f['h2'].value))

        if i == 0:
                n = star_m.shape[1]
                dr = 0.2
                factor = dr*n
                bins = np.arange(0., factor, dr)

                h1_frac_tukey = np.zeros((len(bin_labels), n)); h1_frac_large_scale = np.zeros((len(bin_labels), n)); h1_frac_small_scale = np.zeros((len(bin_labels), n))
                h1_mass_tukey = np.zeros((len(bin_labels), n)); h1_mass_large_scale = np.zeros((len(bin_labels), n)); h1_mass_small_scale = np.zeros((len(bin_labels), n))
                h1_mass_sm_tukey = np.zeros((len(bin_labels), n)); h1_mass_sm_large_scale = np.zeros((len(bin_labels), n)); h1_mass_sm_small_scale = np.zeros((len(bin_labels), n))
            
                h2_frac_tukey = np.zeros((len(bin_labels), n)); h2_frac_large_scale = np.zeros((len(bin_labels), n)); h2_frac_small_scale = np.zeros((len(bin_labels), n))
                h2_mass_tukey = np.zeros((len(bin_labels), n)); h2_mass_large_scale = np.zeros((len(bin_labels), n)); h2_mass_small_scale = np.zeros((len(bin_labels), n))
                h2_mass_sm_tukey = np.zeros((len(bin_labels), n)); h2_mass_sm_large_scale = np.zeros((len(bin_labels), n)); h2_mass_sm_small_scale = np.zeros((len(bin_labels), n))

    	h1_mass_sm = gas_h1*gas_m / star_m
        h2_mass_sm = gas_h2*gas_m / star_m

	no_gals[i] = len(star_m)

        # hi and h2:

        tukey, scale = tukey_biweight(gas_h1)
	h1_frac_tukey[i] = np.log10(tukey)
        h1_frac_large_scale[i] = scale / (np.log(10.)*tukey)
        h1_frac_small_scale[i] = scale / (np.sqrt(no_gals[i])* np.log(10.)*tukey)
        tukey, scale = tukey_biweight(gas_h1*gas_m)
        h1_mass_tukey[i] = np.log10(tukey)
        h1_mass_large_scale[i] = scale / (np.log(10.)*tukey)
        h1_mass_small_scale[i] = scale / (np.sqrt(no_gals[i])* np.log(10.)*tukey)
        tukey, scale = tukey_biweight(h1_mass_sm)
        h1_mass_sm_tukey[i] = np.log10(tukey)
        h1_mass_sm_large_scale[i] = scale / (np.log(10.)*tukey)
        h1_mass_sm_small_scale[i] = scale / (np.sqrt(no_gals[i])* np.log(10.)*tukey)

        tukey, scale = tukey_biweight(gas_h2)
        h2_frac_tukey[i] = np.log10(tukey)
        h2_frac_large_scale[i] = scale / (np.log(10.)*tukey)
        h2_frac_small_scale[i] = scale / (np.sqrt(no_gals[i])* np.log(10.)*tukey)
        tukey, scale = tukey_biweight(gas_h2*gas_m)
        h2_mass_tukey[i] = np.log10(tukey)
        h2_mass_large_scale[i] = scale / (np.log(10.)*tukey)
        h2_mass_small_scale[i] = scale / (np.sqrt(no_gals[i])* np.log(10.)*tukey)
        tukey, scale = tukey_biweight(h2_mass_sm)
        h2_mass_sm_tukey[i] = np.log10(tukey)
        h2_mass_sm_large_scale[i] = scale / (np.log(10.)*tukey)
        h2_mass_sm_small_scale[i] = scale / (np.sqrt(no_gals[i])* np.log(10.)*tukey)


        with h5py.File(results_dir+'h1_frac_data.h5', 'a') as f:
                f.create_dataset('no_gals_'+bin_labels[i], data=np.array(no_gals[i]))
                f.create_dataset('tukey_'+bin_labels[i], data=np.array(h1_frac_tukey[i]))
                f.create_dataset('large_scale_'+bin_labels[i], data=np.array(h1_frac_large_scale[i]))
                f.create_dataset('small_scale_'+bin_labels[i], data=np.array(h1_frac_small_scale[i]))

        with h5py.File(results_dir+'h1_mass_data.h5', 'a') as f:
                f.create_dataset('no_gals_'+bin_labels[i], data=np.array(no_gals[i]))
                f.create_dataset('tukey_'+bin_labels[i], data=np.array(h1_mass_tukey[i]))
                f.create_dataset('large_scale_'+bin_labels[i], data=np.array(h1_mass_large_scale[i]))
                f.create_dataset('small_scale_'+bin_labels[i], data=np.array(h1_mass_small_scale[i]))
                
        with h5py.File(results_dir+'h1_mass_sm_data.h5', 'a') as f:
                f.create_dataset('no_gals_'+bin_labels[i], data=np.array(no_gals[i]))
                f.create_dataset('tukey_'+bin_labels[i], data=np.array(h1_mass_sm_tukey[i]))
                f.create_dataset('large_scale_'+bin_labels[i], data=np.array(h1_mass_sm_large_scale[i]))
                f.create_dataset('small_scale_'+bin_labels[i], data=np.array(h1_mass_sm_small_scale[i]))
                
        with h5py.File(results_dir+'h2_frac_data.h5', 'a') as f:
                f.create_dataset('no_gals_'+bin_labels[i], data=np.array(no_gals[i]))
                f.create_dataset('tukey_'+bin_labels[i], data=np.array(h2_frac_tukey[i]))
                f.create_dataset('large_scale_'+bin_labels[i], data=np.array(h2_frac_large_scale[i]))
                f.create_dataset('small_scale_'+bin_labels[i], data=np.array(h2_frac_small_scale[i]))
                
        with h5py.File(results_dir+'h2_mass_data.h5', 'a') as f:
                f.create_dataset('no_gals_'+bin_labels[i], data=np.array(no_gals[i]))
                f.create_dataset('tukey_'+bin_labels[i], data=np.array(h2_mass_tukey[i]))
                f.create_dataset('large_scale_'+bin_labels[i], data=np.array(h2_mass_large_scale[i]))
                f.create_dataset('small_scale_'+bin_labels[i], data=np.array(h2_mass_small_scale[i]))
                
        with h5py.File(results_dir+'h2_mass_sm_data.h5', 'a') as f:
                f.create_dataset('no_gals_'+bin_labels[i], data=np.array(no_gals[i]))
                f.create_dataset('tukey_'+bin_labels[i], data=np.array(h2_mass_sm_tukey[i]))
                f.create_dataset('large_scale_'+bin_labels[i], data=np.array(h2_mass_sm_large_scale[i]))
                f.create_dataset('small_scale_'+bin_labels[i], data=np.array(h2_mass_sm_small_scale[i]))
