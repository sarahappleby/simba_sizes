from plotting_methods import tukey_biweight
import h5py
import numpy as np
import matplotlib.pyplot as plt

model = 'm100n1024'
wind = 's50j7k'
snap = '151'
selection = 'sf'
m = 3
bin_label = '10.5-11.0'

epsilon = [1.e-4, 1.e-5, 1.e-7, 1.e-9, 1.e-11, 1.e-13, 1.e-15]

if selection =='gv':
    name = 'green_valley'
elif selection == 'sf':
    name = 'star_forming'
basic_dir = '/home/sapple/simba_sizes/profiles/paper/'

centrals_dir = basic_dir + 'centrals/'+model+'_'+snap+'/'+wind+'/'+name+'/random_orientation/'
results_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/data/'

with h5py.File(centrals_dir+'mask_'+str(m)+'_all_profiles.h5', 'r') as f:
        cen_star_m = f['sm'].value
        cen_gas_sfr = f['gas_sfr'].value
        cen_gas_h1 = f['h1_m'].value
        cen_gas_h2 = f['h2_m'].value

n = cen_star_m.shape[1]
cen_no_gals = cen_star_m.shape[0]
cen_ssfr_tukey = np.zeros((len(epsilon), n)); cen_ssfr_large_scale = np.zeros((len(epsilon), n)); cen_ssfr_small_scale = np.zeros((len(epsilon), n))
cen_sfr_tukey = np.zeros((len(epsilon), n)); cen_sfr_large_scale = np.zeros((len(epsilon), n)); cen_sfr_small_scale = np.zeros((len(epsilon), n))
cen_h1_tukey = np.zeros((len(epsilon), n)); cen_h1_large_scale = np.zeros((len(epsilon), n)); cen_h1_small_scale = np.zeros((len(epsilon), n))
cen_fmol_tukey = np.zeros((len(epsilon), n)); cen_fmol_large_scale = np.zeros((len(epsilon), n)); cen_fmol_small_scale = np.zeros((len(epsilon), n))

cen_gas_ssfr = cen_gas_sfr / cen_star_m
cen_gas_fmol = cen_gas_h2 / cen_gas_h1

dr = 0.2
factor = dr*n
bins = np.arange(0., factor, dr)
rplot = bins+(dr*0.5)

fig, ax = plt.subplots(2, 2, figsize=(13, 13))
ax = ax.flatten()

med_sfr = np.nanmedian(cen_gas_sfr, axis=0)
med_sfr[np.where(med_sfr == 0.)[0]] = 1.e-6
med_sfr = np.log10(med_sfr)
med_ssfr = np.nanmedian(cen_gas_ssfr, axis=0)
med_ssfr[np.where(med_ssfr == 0.)] = 1.e-14
med_ssfr = np.log10(med_ssfr)
med_h1 = np.nanmedian(np.log10(cen_gas_h1), axis=0)
med_fmol = np.nanmedian(np.log10(cen_gas_fmol), axis=0)

ax[0].plot(rplot, med_sfr, linestyle='--', label='median')
ax[1].plot(rplot, med_ssfr, linestyle='--', color='k')
ax[2].plot(rplot, med_h1, linestyle='--', )
ax[3].plot(rplot, med_fmol, linestyle='--', )


for i, e in enumerate(epsilon):

    tukey, scale = tukey_biweight(cen_gas_sfr, epsilon=e)
    tukey[np.where(tukey == 0.)[0]] = 1.e-6
    cen_sfr_tukey[i] = np.log10(tukey)
    cen_sfr_large_scale[i] = scale / (np.log(10.)*tukey)
    cen_sfr_small_scale[i] = scale / (np.sqrt(cen_no_gals)* np.log(10.)*tukey)

    tukey, scale = tukey_biweight(cen_gas_ssfr, epsilon=e)
    cen_ssfr_tukey[i] = np.log10(tukey)
    cen_ssfr_large_scale[i] = scale / (np.log(10.)*tukey)
    cen_ssfr_small_scale[i] = scale / (np.sqrt(cen_no_gals)* np.log(10.)*tukey)

    tukey, scale = tukey_biweight(cen_gas_h1, epsilon=e)
    cen_h1_tukey[i] = np.log10(tukey)
    cen_h1_large_scale[i] = scale / (np.log(10.)*tukey)
    cen_h1_small_scale[i] = scale / (np.sqrt(cen_no_gals)* np.log(10.)*tukey)

    tukey, scale = tukey_biweight(cen_gas_fmol, epsilon=e)
    tukey[np.where(tukey == 0.)[0]] = 1.e-6
    cen_fmol_tukey[i] = np.log10(tukey)
    cen_fmol_large_scale[i] = scale / (np.log(10.)*tukey)
    cen_fmol_small_scale[i] = scale / (np.sqrt(cen_no_gals)* np.log(10.)*tukey)


    ax[0].plot(rplot, cen_sfr_tukey[i], linestyle='--', label='epsilon = '+str(e)) 
    ax[1].plot(rplot, cen_ssfr_tukey[i], linestyle='--', label='epsilon = '+str(e))
    ax[2].plot(rplot, cen_h1_tukey[i], linestyle='--', label='epsilon = '+str(e))
    ax[3].plot(rplot, cen_fmol_tukey[i], linestyle='--', label='epsilon = '+str(e))

ax[0].set_ylabel('SFR')
ax[1].set_ylabel('sSFR')
ax[2].set_ylabel('HI')
ax[3].set_ylabel('fmol')
plt.legend()
plt.savefig(results_dir+'epsilon_test.png')
plt.show()
plt.clf()


