import h5py
import numpy as np
import sys
import caesar
from plotting_methods import *

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

mass_bins = [10., 10.5, 11.0]
bin_labels = ['10.0-10.5', '10.5-11.0', '>11.0']
#mass_bins = [10., 10.6]
#bin_labels = ['10.0-10.6','>10.6']


basic_dir = '/home/sapple/simba_sizes/profiles/paper/with_dust_sizes/'

profs_file = basic_dir + '/'+wind+'/'+model+'_'+snap+'/all_profiles_'+angle+'.h5'
#profs_file = basic_dir + '/'+wind+'/resolution_test/'+model+'_'+snap+'/all_profiles_'+angle+'.h5'
results_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/data/'
if not wind == 's50':
    results_dir += wind+'_'
if model == 'm25n512':
    results_dir += model + '_'

if angle == 'rotated_faceon':
    angle = 'rot'
elif angle == 'random_orientation':
    angle = 'rand'

results_dir += selection + '_' + angle

data_dir = '/home/rad/data/'+model+'/'+wind+'/'
sim =  caesar.load(data_dir+'Groups/'+model+'_'+snap+'.hdf5', LoadHalo=False)
gal_sm = np.log10(np.array([i.masses['stellar'].in_units('Msun') for i in sim.galaxies]))
gal_cent = np.array([i.central for i in sim.galaxies])
gal_ids = np.array([False for i in range(len(sim.galaxies))])
gal_ids[gals] = True

# SA 19/2/20: removing sample_mask - it's identical to gals??
#with h5py.File(sample_file, 'r') as sf:
#   sample_mask = sf[model+'_'+snap][:]


cen_no_gals = np.zeros(len(bin_labels))
sat_no_gals = np.zeros(len(bin_labels))
all_no_gals = np.zeros(len(bin_labels))

for i, b in enumerate(bin_labels):

    if i != len(mass_bins) -1:
        mass_mask = (gal_sm > mass_bins[i]) & (gal_sm < mass_bins[i+1])
    else:
        mass_mask = gal_sm > mass_bins[i]

    #mask = mass_mask*gal_ids*gal_cent*sample_mask
    mask = mass_mask*gal_ids*gal_cent
    with h5py.File(profs_file, 'r') as f:
        cen_star_m = f['sm'][:][mask]
        cen_gas_sfr = f['gas_sfr'][:][mask]
        cen_gas_h1 = f['h1_m'][:][mask]
        cen_gas_h2 = f['h2_m'][:][mask]
        cen_ngas = f['gas_npart'][:][mask]
        cen_nstar = f['star_npart'][:][mask]

    #mask = mass_mask*gal_ids*np.invert(gal_cent)*sample_mask
    mask = mass_mask*gal_ids*np.invert(gal_cent)
    with h5py.File(profs_file, 'r') as f:
        sat_star_m = f['sm'][:][mask]
        sat_gas_sfr = f['gas_sfr'][:][mask]
        sat_gas_h1 = f['h1_m'][:][mask]
        sat_gas_h2 = f['h2_m'][:][mask] 
        sat_ngas = f['gas_npart'][:][mask]
        sat_nstar = f['star_npart'][:][mask]
                
                
    print(str(len(sat_star_m)), str(len(cen_star_m)))

    if i == 0:
        n = cen_star_m.shape[1]
        cen_ssfr_tukey = np.zeros((len(bin_labels), n)); cen_ssfr_large_scale = np.zeros((len(bin_labels), n)); cen_ssfr_small_scale = np.zeros((len(bin_labels), n))
        cen_sfr_tukey = np.zeros((len(bin_labels), n)); cen_sfr_large_scale = np.zeros((len(bin_labels), n)); cen_sfr_small_scale = np.zeros((len(bin_labels), n))
        cen_h1_tukey = np.zeros((len(bin_labels), n)); cen_h1_large_scale = np.zeros((len(bin_labels), n)); cen_h1_small_scale = np.zeros((len(bin_labels), n))
        cen_h2_tukey = np.zeros((len(bin_labels), n)); cen_h2_large_scale = np.zeros((len(bin_labels), n)); cen_h2_small_scale = np.zeros((len(bin_labels), n))
        cen_fmol_tukey = np.zeros((len(bin_labels), n)); cen_fmol_large_scale = np.zeros((len(bin_labels), n)); cen_fmol_small_scale = np.zeros((len(bin_labels), n))
        cen_sfe_tukey = np.zeros((len(bin_labels), n)); cen_sfe_large_scale = np.zeros((len(bin_labels), n)); cen_sfe_small_scale = np.zeros((len(bin_labels), n))
        cen_fh2_tukey = np.zeros((len(bin_labels), n)); cen_fh2_large_scale = np.zeros((len(bin_labels), n)); cen_fh2_small_scale = np.zeros((len(bin_labels), n))
        cen_ngas_tukey = np.zeros((len(bin_labels), n)); cen_ngas_large_scale = np.zeros((len(bin_labels), n)); cen_ngas_small_scale = np.zeros((len(bin_labels), n))
        cen_nstar_tukey = np.zeros((len(bin_labels), n)); cen_nstar_large_scale = np.zeros((len(bin_labels), n)); cen_nstar_small_scale = np.zeros((len(bin_labels), n))

        sat_ssfr_tukey = np.zeros((len(bin_labels), n)); sat_ssfr_large_scale = np.zeros((len(bin_labels), n)); sat_ssfr_small_scale = np.zeros((len(bin_labels), n))
        sat_sfr_tukey = np.zeros((len(bin_labels), n)); sat_sfr_large_scale = np.zeros((len(bin_labels), n)); sat_sfr_small_scale = np.zeros((len(bin_labels), n))
        sat_h1_tukey = np.zeros((len(bin_labels), n)); sat_h1_large_scale = np.zeros((len(bin_labels), n)); sat_h1_small_scale = np.zeros((len(bin_labels), n))
        sat_h2_tukey = np.zeros((len(bin_labels), n)); sat_h2_large_scale = np.zeros((len(bin_labels), n)); sat_h2_small_scale = np.zeros((len(bin_labels), n))
        sat_fmol_tukey = np.zeros((len(bin_labels), n)); sat_fmol_large_scale = np.zeros((len(bin_labels), n)); sat_fmol_small_scale = np.zeros((len(bin_labels), n))
        sat_sfe_tukey = np.zeros((len(bin_labels), n)); sat_sfe_large_scale = np.zeros((len(bin_labels), n)); sat_sfe_small_scale = np.zeros((len(bin_labels), n))
        sat_fh2_tukey = np.zeros((len(bin_labels), n)); sat_fh2_large_scale = np.zeros((len(bin_labels), n)); sat_fh2_small_scale = np.zeros((len(bin_labels), n))
        sat_ngas_tukey = np.zeros((len(bin_labels), n)); sat_ngas_large_scale = np.zeros((len(bin_labels), n)); sat_ngas_small_scale = np.zeros((len(bin_labels), n))
        sat_nstar_tukey = np.zeros((len(bin_labels), n)); sat_nstar_large_scale = np.zeros((len(bin_labels), n)); sat_nstar_small_scale = np.zeros((len(bin_labels), n))

        all_ssfr_tukey = np.zeros((len(bin_labels), n)); all_ssfr_large_scale = np.zeros((len(bin_labels), n)); all_ssfr_small_scale = np.zeros((len(bin_labels), n))
        all_sfr_tukey = np.zeros((len(bin_labels), n)); all_sfr_large_scale = np.zeros((len(bin_labels), n)); all_sfr_small_scale = np.zeros((len(bin_labels), n))
        all_h1_tukey = np.zeros((len(bin_labels), n)); all_h1_large_scale = np.zeros((len(bin_labels), n)); all_h1_small_scale = np.zeros((len(bin_labels), n))
        all_h2_tukey = np.zeros((len(bin_labels), n)); all_h2_large_scale = np.zeros((len(bin_labels), n)); all_h2_small_scale = np.zeros((len(bin_labels), n))
        all_fmol_tukey = np.zeros((len(bin_labels), n)); all_fmol_large_scale = np.zeros((len(bin_labels), n)); all_fmol_small_scale = np.zeros((len(bin_labels), n))
        all_sfe_tukey = np.zeros((len(bin_labels), n)); all_sfe_large_scale = np.zeros((len(bin_labels), n)); all_sfe_small_scale = np.zeros((len(bin_labels), n))
        all_fh2_tukey = np.zeros((len(bin_labels), n)); all_fh2_large_scale = np.zeros((len(bin_labels), n)); all_fh2_small_scale = np.zeros((len(bin_labels), n))
        all_ngas_tukey = np.zeros((len(bin_labels), n)); all_ngas_large_scale = np.zeros((len(bin_labels), n)); all_ngas_small_scale = np.zeros((len(bin_labels), n))
        all_nstar_tukey = np.zeros((len(bin_labels), n)); all_nstar_large_scale = np.zeros((len(bin_labels), n)); all_nstar_small_scale = np.zeros((len(bin_labels), n))
        
    # centrals:
    cen_no_gals[i] = len(cen_gas_sfr)
    cen_gas_sfe = cen_gas_sfr / cen_gas_h2
    cen_gas_fh2 = cen_gas_h2 / cen_star_m
    cen_gas_ssfr = cen_gas_sfr / cen_star_m
    cen_gas_fmol = cen_gas_h2 / cen_gas_h1

    tukey, scale = tukey_biweight(cen_gas_sfr)
    tukey[np.where(tukey == 0.)[0]] = 1.e-6
    cen_sfr_tukey[i] = np.log10(tukey)
    cen_sfr_large_scale[i] = scale / (np.log(10.)*tukey) 
    cen_sfr_small_scale[i] = scale / (np.sqrt(cen_no_gals[i])* np.log(10.)*tukey)

    tukey, scale = tukey_biweight(cen_gas_ssfr)
    tukey[np.where(tukey == 0.)[0]] = 10.**-13.5
    cen_ssfr_tukey[i] = np.log10(tukey)
    cen_ssfr_large_scale[i] = scale / (np.log(10.)*tukey)
    cen_ssfr_small_scale[i] = scale / (np.sqrt(cen_no_gals[i])* np.log(10.)*tukey)

    tukey, scale = tukey_biweight(cen_gas_h1)
    cen_h1_tukey[i] = np.log10(tukey)
    cen_h1_large_scale[i] = scale / (np.log(10.)*tukey)
    cen_h1_small_scale[i] = scale / (np.sqrt(cen_no_gals[i])* np.log(10.)*tukey)

    tukey, scale = tukey_biweight(cen_gas_h2)
    cen_h2_tukey[i] = np.log10(tukey)
    cen_h2_large_scale[i] = scale / (np.log(10.)*tukey)
    cen_h2_small_scale[i] = scale / (np.sqrt(cen_no_gals[i])* np.log(10.)*tukey)

    tukey, scale = tukey_biweight(cen_gas_fmol)
    tukey[np.where(tukey == 0.)[0]] = 1.e-6
    if snap == '151':
            for j in range(len(scale) -1):
                if scale[j+1] > 2.*scale[j]:
                    scale[j+1] = np.median(scale)
    cen_fmol_tukey[i] = np.log10(tukey)
    cen_fmol_large_scale[i] = scale / (np.log(10.)*tukey)
    cen_fmol_small_scale[i] = scale / (np.sqrt(cen_no_gals[i])* np.log(10.)*tukey)

    tukey, scale = tukey_biweight(cen_gas_sfe)
    cen_sfe_tukey[i] = np.log10(tukey)
    cen_sfe_large_scale[i] = scale / (np.log(10.)*tukey)
    cen_sfe_small_scale[i] = scale / (np.sqrt(cen_no_gals[i])* np.log(10.)*tukey)

    tukey, scale = tukey_biweight(cen_gas_fh2)
    cen_fh2_tukey[i] = np.log10(tukey)
    cen_fh2_large_scale[i] = scale / (np.log(10.)*tukey)
    cen_fh2_small_scale[i] = scale / (np.sqrt(cen_no_gals[i])* np.log(10.)*tukey)

    tukey, scale = tukey_biweight(cen_ngas)
    cen_ngas_tukey[i] = np.log10(tukey)
    cen_ngas_large_scale[i] = scale / (np.log(10.)*tukey)
    cen_ngas_small_scale[i] = scale / (np.sqrt(cen_no_gals[i])* np.log(10.)*tukey)

    tukey, scale = tukey_biweight(cen_nstar)
    cen_nstar_tukey[i] = np.log10(tukey)
    cen_nstar_large_scale[i] = scale / (np.log(10.)*tukey)
    cen_nstar_small_scale[i] = scale / (np.sqrt(cen_no_gals[i])* np.log(10.)*tukey)

    # satellites:
    sat_no_gals[i] = len(sat_gas_sfr)
    sat_gas_sfe = sat_gas_sfr / sat_gas_h2
    sat_gas_fh2 = sat_gas_h2 / sat_star_m
    sat_gas_ssfr = sat_gas_sfr / sat_star_m
    sat_gas_fmol = sat_gas_h2 / sat_gas_h1

    tukey, scale = tukey_biweight(sat_gas_sfr)
    tukey[np.where(tukey == 0.)[0]] = 1.e-6
    sat_sfr_tukey[i] = np.log10(tukey)
    sat_sfr_large_scale[i] = scale / (np.log(10.)*tukey)
    sat_sfr_small_scale[i] = scale / (np.sqrt(sat_no_gals[i])* np.log(10.)*tukey)

    tukey, scale = tukey_biweight(sat_gas_ssfr)
    tukey[np.where(tukey == 0.)[0]] = 10.**-13.5
    sat_ssfr_tukey[i] = np.log10(tukey)
    sat_ssfr_large_scale[i] = scale / (np.log(10.)*tukey)
    sat_ssfr_small_scale[i] = scale / (np.sqrt(sat_no_gals[i])* np.log(10.)*tukey)

    tukey, scale = tukey_biweight(sat_gas_h1)
    sat_h1_tukey[i] = np.log10(tukey)
    sat_h1_large_scale[i] = scale / (np.log(10.)*tukey)
    sat_h1_small_scale[i] = scale / (np.sqrt(cen_no_gals[i])* np.log(10.)*tukey)

    tukey, scale = tukey_biweight(sat_gas_h2)
    sat_h2_tukey[i] = np.log10(tukey)
    sat_h2_large_scale[i] = scale / (np.log(10.)*tukey)
    sat_h2_small_scale[i] = scale / (np.sqrt(cen_no_gals[i])* np.log(10.)*tukey)

    tukey, scale = tukey_biweight(sat_gas_fmol)
    tukey[np.where(tukey == 0.)[0]] = 1.e-6
    sat_fmol_tukey[i] = np.log10(tukey)
    sat_fmol_large_scale[i] = scale / (np.log(10.)*tukey)
    sat_fmol_small_scale[i] = scale / (np.sqrt(sat_no_gals[i])* np.log(10.)*tukey)

    tukey, scale = tukey_biweight(sat_gas_sfe)
    sat_sfe_tukey[i] = np.log10(tukey)
    sat_sfe_large_scale[i] = scale / (np.log(10.)*tukey)
    sat_sfe_small_scale[i] = scale / (np.sqrt(sat_no_gals[i])* np.log(10.)*tukey)

    tukey, scale = tukey_biweight(sat_gas_fh2)
    sat_fh2_tukey[i] = np.log10(tukey)
    sat_fh2_large_scale[i] = scale / (np.log(10.)*tukey)
    sat_fh2_small_scale[i] = scale / (np.sqrt(sat_no_gals[i])* np.log(10.)*tukey)

    tukey, scale = tukey_biweight(sat_ngas)
    sat_ngas_tukey[i] = np.log10(tukey)
    sat_ngas_large_scale[i] = scale / (np.log(10.)*tukey)
    sat_ngas_small_scale[i] = scale / (np.sqrt(sat_no_gals[i])* np.log(10.)*tukey)

    tukey, scale = tukey_biweight(sat_nstar)
    sat_nstar_tukey[i] = np.log10(tukey)
    sat_nstar_large_scale[i] = scale / (np.log(10.)*tukey)
    sat_nstar_small_scale[i] = scale / (np.sqrt(sat_no_gals[i])* np.log(10.)*tukey)

    # all:
    ssfr = np.concatenate((cen_gas_ssfr, sat_gas_ssfr))
    sfr = np.concatenate((cen_gas_sfr, sat_gas_sfr))
    h1 = np.concatenate((cen_gas_h1, sat_gas_h1))
    h2 = np.concatenate((cen_gas_h2, sat_gas_h2))
    fmol = np.concatenate((cen_gas_fmol, sat_gas_fmol))
    sfe = np.concatenate((cen_gas_sfe, sat_gas_sfe))
    fh2 = np.concatenate((cen_gas_fh2, sat_gas_fh2))
    ngas = np.concatenate((cen_ngas, sat_ngas))
    nstar = np.concatenate((cen_nstar, sat_nstar))

    all_no_gals[i] = len(sfr)

    tukey, scale = tukey_biweight(ssfr)
    all_ssfr_tukey[i] = np.log10(tukey)
    all_ssfr_large_scale[i] = scale / (np.log(10.)*tukey)
    all_ssfr_small_scale[i] = scale / (np.sqrt(all_no_gals[i])* np.log(10.)*tukey)

    tukey, scale = tukey_biweight(sfr)
    tukey[np.where(tukey == 0.)[0]] = 1.e-6
    all_sfr_tukey[i] = np.log10(tukey)
    all_sfr_large_scale[i] = scale / (np.log(10.)*tukey)
    all_sfr_small_scale[i] = scale / (np.sqrt(all_no_gals[i])* np.log(10.)*tukey)

    tukey, scale = tukey_biweight(h1)
    tukey[np.where(tukey == 0.)[0]] = 10.**-13.5
    all_h1_tukey[i] = np.log10(tukey)
    all_h1_large_scale[i] = scale / (np.log(10.)*tukey)
    all_h1_small_scale[i] = scale / (np.sqrt(all_no_gals[i])* np.log(10.)*tukey)

    tukey, scale = tukey_biweight(h2)
    all_h2_tukey[i] = np.log10(tukey)
    all_h2_large_scale[i] = scale / (np.log(10.)*tukey)
    all_h2_small_scale[i] = scale / (np.sqrt(all_no_gals[i])* np.log(10.)*tukey)

    tukey, scale = tukey_biweight(fmol)
    tukey[np.where(tukey == 0.)[0]] = 1.e-6
    if snap == '151':
        for j in range(len(scale) -1):
            if scale[j+1] > 2.*scale[j]:
                scale[j+1] = np.median(scale)
    all_fmol_tukey[i] = np.log10(tukey)
    all_fmol_large_scale[i] = scale / (np.log(10.)*tukey)
    all_fmol_small_scale[i] = scale / (np.sqrt(all_no_gals[i])* np.log(10.)*tukey)

    tukey, scale = tukey_biweight(sfe)
    all_sfe_tukey[i] = np.log10(tukey)
    all_sfe_large_scale[i] = scale / (np.log(10.)*tukey)
    all_sfe_small_scale[i] = scale / (np.sqrt(all_no_gals[i])* np.log(10.)*tukey)

    tukey, scale = tukey_biweight(fh2)
    tukey[np.where(tukey == 0.)[0]] = 1.e-6
    all_fh2_tukey[i] = np.log10(tukey)
    all_fh2_large_scale[i] = scale / (np.log(10.)*tukey)
    all_fh2_small_scale[i] = scale / (np.sqrt(all_no_gals[i])* np.log(10.)*tukey)
    
    tukey, scale = tukey_biweight(ngas)
    all_ngas_tukey[i] = np.log10(tukey)
    all_ngas_large_scale[i] = scale / (np.log(10.)*tukey)
    all_ngas_small_scale[i] = scale / (np.sqrt(all_no_gals[i])* np.log(10.)*tukey)

    tukey, scale = tukey_biweight(nstar)
    all_nstar_tukey[i] = np.log10(tukey)
    all_nstar_large_scale[i] = scale / (np.log(10.)*tukey)
    all_nstar_small_scale[i] = scale / (np.sqrt(all_no_gals[i])* np.log(10.)*tukey)

    """
    with h5py.File(results_dir+'_ssfr_data.h5', 'a') as f:
        f.create_dataset('cen_no_gals_'+bin_labels[i], data=np.array(cen_no_gals[i]))
        f.create_dataset('cen_tukey_'+bin_labels[i], data=np.array(cen_ssfr_tukey[i]))
        f.create_dataset('cen_large_scale_'+bin_labels[i], data=np.array(cen_ssfr_large_scale[i]))
        f.create_dataset('cen_small_scale_'+bin_labels[i], data=np.array(cen_ssfr_small_scale[i]))
        f.create_dataset('all_no_gals_'+bin_labels[i], data=np.array(all_no_gals[i]))
        f.create_dataset('sat_no_gals_'+bin_labels[i], data=np.array(sat_no_gals[i]))
        f.create_dataset('all_tukey_'+bin_labels[i], data=np.array(all_ssfr_tukey[i]))
        f.create_dataset('sat_tukey_'+bin_labels[i], data=np.array(sat_ssfr_tukey[i]))
        f.create_dataset('all_large_scale_'+bin_labels[i], data=np.array(all_ssfr_large_scale[i]))
        f.create_dataset('sat_large_scale_'+bin_labels[i], data=np.array(sat_ssfr_large_scale[i]))
        f.create_dataset('all_small_scale_'+bin_labels[i], data=np.array(all_ssfr_small_scale[i]))
        f.create_dataset('sat_small_scale_'+bin_labels[i], data=np.array(sat_ssfr_small_scale[i]))

    with h5py.File(results_dir+'_sfr_data.h5', 'a') as f:
        f.create_dataset('cen_no_gals_'+bin_labels[i], data=np.array(cen_no_gals[i]))
        f.create_dataset('cen_tukey_'+bin_labels[i], data=np.array(cen_sfr_tukey[i]))
        f.create_dataset('cen_large_scale_'+bin_labels[i], data=np.array(cen_sfr_large_scale[i]))
        f.create_dataset('cen_small_scale_'+bin_labels[i], data=np.array(cen_sfr_small_scale[i]))
        f.create_dataset('all_no_gals_'+bin_labels[i], data=np.array(all_no_gals[i]))
        f.create_dataset('all_tukey_'+bin_labels[i], data=np.array(all_sfr_tukey[i]))
        f.create_dataset('all_large_scale_'+bin_labels[i], data=np.array(all_sfr_large_scale[i]))
        f.create_dataset('all_small_scale_'+bin_labels[i], data=np.array(all_sfr_small_scale[i]))
        f.create_dataset('sat_no_gals_'+bin_labels[i], data=np.array(sat_no_gals[i]))
        f.create_dataset('sat_tukey_'+bin_labels[i], data=np.array(sat_sfr_tukey[i]))
        f.create_dataset('sat_large_scale_'+bin_labels[i], data=np.array(sat_sfr_large_scale[i]))
        f.create_dataset('sat_small_scale_'+bin_labels[i], data=np.array(sat_sfr_small_scale[i]))

    with h5py.File(results_dir+'_h1_data.h5', 'a') as f:
        f.create_dataset('cen_no_gals_'+bin_labels[i], data=np.array(cen_no_gals[i]))
        f.create_dataset('cen_tukey_'+bin_labels[i], data=np.array(cen_h1_tukey[i]))
        f.create_dataset('cen_large_scale_'+bin_labels[i], data=np.array(cen_h1_large_scale[i]))
        f.create_dataset('cen_small_scale_'+bin_labels[i], data=np.array(cen_h1_small_scale[i]))
        f.create_dataset('all_no_gals_'+bin_labels[i], data=np.array(all_no_gals[i]))
        f.create_dataset('all_tukey_'+bin_labels[i], data=np.array(all_h1_tukey[i]))
        f.create_dataset('all_large_scale_'+bin_labels[i], data=np.array(all_h1_large_scale[i]))
        f.create_dataset('all_small_scale_'+bin_labels[i], data=np.array(all_h1_small_scale[i]))
        f.create_dataset('sat_no_gals_'+bin_labels[i], data=np.array(sat_no_gals[i]))
        f.create_dataset('sat_tukey_'+bin_labels[i], data=np.array(sat_h1_tukey[i]))
        f.create_dataset('sat_large_scale_'+bin_labels[i], data=np.array(sat_h1_large_scale[i]))
        f.create_dataset('sat_small_scale_'+bin_labels[i], data=np.array(sat_h1_small_scale[i]))
    
    with h5py.File(results_dir+'_h2_data.h5', 'a') as f:
        f.create_dataset('cen_no_gals_'+bin_labels[i], data=np.array(cen_no_gals[i]))
        f.create_dataset('cen_tukey_'+bin_labels[i], data=np.array(cen_h2_tukey[i]))
        f.create_dataset('cen_large_scale_'+bin_labels[i], data=np.array(cen_h2_large_scale[i]))
        f.create_dataset('cen_small_scale_'+bin_labels[i], data=np.array(cen_h2_small_scale[i]))
        f.create_dataset('all_no_gals_'+bin_labels[i], data=np.array(all_no_gals[i]))
        f.create_dataset('all_tukey_'+bin_labels[i], data=np.array(all_h2_tukey[i]))
        f.create_dataset('all_large_scale_'+bin_labels[i], data=np.array(all_h2_large_scale[i]))
        f.create_dataset('all_small_scale_'+bin_labels[i], data=np.array(all_h2_small_scale[i]))
        f.create_dataset('sat_no_gals_'+bin_labels[i], data=np.array(sat_no_gals[i]))
        f.create_dataset('sat_tukey_'+bin_labels[i], data=np.array(sat_h2_tukey[i]))
        f.create_dataset('sat_large_scale_'+bin_labels[i], data=np.array(sat_h2_large_scale[i]))
        f.create_dataset('sat_small_scale_'+bin_labels[i], data=np.array(sat_h2_small_scale[i]))
    
    with h5py.File(results_dir+'_fmol_data.h5', 'a') as f:
        f.create_dataset('cen_no_gals_'+bin_labels[i], data=np.array(cen_no_gals[i]))
        f.create_dataset('cen_tukey_'+bin_labels[i], data=np.array(cen_fmol_tukey[i]))
        f.create_dataset('cen_large_scale_'+bin_labels[i], data=np.array(cen_fmol_large_scale[i]))
        f.create_dataset('cen_small_scale_'+bin_labels[i], data=np.array(cen_fmol_small_scale[i]))
        f.create_dataset('all_no_gals_'+bin_labels[i], data=np.array(all_no_gals[i]))
        f.create_dataset('all_tukey_'+bin_labels[i], data=np.array(all_fmol_tukey[i]))
        f.create_dataset('all_large_scale_'+bin_labels[i], data=np.array(all_fmol_large_scale[i]))
        f.create_dataset('all_small_scale_'+bin_labels[i], data=np.array(all_fmol_small_scale[i]))
        f.create_dataset('sat_no_gals_'+bin_labels[i], data=np.array(sat_no_gals[i]))
        f.create_dataset('sat_tukey_'+bin_labels[i], data=np.array(sat_fmol_tukey[i]))
        f.create_dataset('sat_large_scale_'+bin_labels[i], data=np.array(sat_fmol_large_scale[i]))
        f.create_dataset('sat_small_scale_'+bin_labels[i], data=np.array(sat_fmol_small_scale[i]))
    
    with h5py.File(results_dir+'_sfe_data.h5', 'a') as f:
        f.create_dataset('cen_no_gals_'+bin_labels[i], data=np.array(cen_no_gals[i]))
        f.create_dataset('cen_tukey_'+bin_labels[i], data=np.array(cen_sfe_tukey[i]))
        f.create_dataset('cen_large_scale_'+bin_labels[i], data=np.array(cen_sfe_large_scale[i]))
        f.create_dataset('cen_small_scale_'+bin_labels[i], data=np.array(cen_sfe_small_scale[i]))
        f.create_dataset('all_no_gals_'+bin_labels[i], data=np.array(all_no_gals[i]))
        f.create_dataset('sat_no_gals_'+bin_labels[i], data=np.array(sat_no_gals[i]))
        f.create_dataset('all_tukey_'+bin_labels[i], data=np.array(all_sfe_tukey[i]))
        f.create_dataset('sat_tukey_'+bin_labels[i], data=np.array(sat_sfe_tukey[i]))
        f.create_dataset('all_large_scale_'+bin_labels[i], data=np.array(all_sfe_large_scale[i]))
        f.create_dataset('sat_large_scale_'+bin_labels[i], data=np.array(sat_sfe_large_scale[i]))
        f.create_dataset('all_small_scale_'+bin_labels[i], data=np.array(all_sfe_small_scale[i]))
        f.create_dataset('sat_small_scale_'+bin_labels[i], data=np.array(sat_sfe_small_scale[i]))

    with h5py.File(results_dir+'_fh2_data.h5', 'a') as f:
        f.create_dataset('cen_no_gals_'+bin_labels[i], data=np.array(cen_no_gals[i]))
        f.create_dataset('cen_tukey_'+bin_labels[i], data=np.array(cen_fh2_tukey[i]))
        f.create_dataset('cen_large_scale_'+bin_labels[i], data=np.array(cen_fh2_large_scale[i]))
        f.create_dataset('cen_small_scale_'+bin_labels[i], data=np.array(cen_fh2_small_scale[i]))
        f.create_dataset('all_no_gals_'+bin_labels[i], data=np.array(all_no_gals[i]))
        f.create_dataset('all_tukey_'+bin_labels[i], data=np.array(all_fh2_tukey[i]))
        f.create_dataset('all_large_scale_'+bin_labels[i], data=np.array(all_fh2_large_scale[i]))
        f.create_dataset('all_small_scale_'+bin_labels[i], data=np.array(all_fh2_small_scale[i]))
        f.create_dataset('sat_no_gals_'+bin_labels[i], data=np.array(sat_no_gals[i]))
        f.create_dataset('sat_tukey_'+bin_labels[i], data=np.array(sat_fh2_tukey[i]))
        f.create_dataset('sat_large_scale_'+bin_labels[i], data=np.array(sat_fh2_large_scale[i]))
        f.create_dataset('sat_small_scale_'+bin_labels[i], data=np.array(sat_fh2_small_scale[i]))
    """

    with h5py.File(results_dir+'_ngas_data.h5', 'a') as f:
            f.create_dataset('cen_no_gals_'+bin_labels[i], data=np.array(cen_no_gals[i]))
            f.create_dataset('cen_tukey_'+bin_labels[i], data=np.array(cen_ngas_tukey[i]))
            f.create_dataset('cen_large_scale_'+bin_labels[i], data=np.array(cen_ngas_large_scale[i]))
            f.create_dataset('cen_small_scale_'+bin_labels[i], data=np.array(cen_ngas_small_scale[i]))
            f.create_dataset('all_no_gals_'+bin_labels[i], data=np.array(all_no_gals[i]))
            f.create_dataset('all_tukey_'+bin_labels[i], data=np.array(all_ngas_tukey[i]))
            f.create_dataset('all_large_scale_'+bin_labels[i], data=np.array(all_ngas_large_scale[i]))
            f.create_dataset('all_small_scale_'+bin_labels[i], data=np.array(all_ngas_small_scale[i]))
            f.create_dataset('sat_no_gals_'+bin_labels[i], data=np.array(sat_no_gals[i]))
            f.create_dataset('sat_tukey_'+bin_labels[i], data=np.array(sat_ngas_tukey[i]))
            f.create_dataset('sat_large_scale_'+bin_labels[i], data=np.array(sat_ngas_large_scale[i]))
            f.create_dataset('sat_small_scale_'+bin_labels[i], data=np.array(sat_ngas_small_scale[i]))

    with h5py.File(results_dir+'_nstar_data.h5', 'a') as f:
            f.create_dataset('cen_no_gals_'+bin_labels[i], data=np.array(cen_no_gals[i]))
            f.create_dataset('cen_tukey_'+bin_labels[i], data=np.array(cen_nstar_tukey[i]))
            f.create_dataset('cen_large_scale_'+bin_labels[i], data=np.array(cen_nstar_large_scale[i]))
            f.create_dataset('cen_small_scale_'+bin_labels[i], data=np.array(cen_nstar_small_scale[i]))
            f.create_dataset('all_no_gals_'+bin_labels[i], data=np.array(all_no_gals[i]))
            f.create_dataset('all_tukey_'+bin_labels[i], data=np.array(all_nstar_tukey[i]))
            f.create_dataset('all_large_scale_'+bin_labels[i], data=np.array(all_nstar_large_scale[i]))
            f.create_dataset('all_small_scale_'+bin_labels[i], data=np.array(all_nstar_small_scale[i]))
            f.create_dataset('sat_no_gals_'+bin_labels[i], data=np.array(sat_no_gals[i]))
            f.create_dataset('sat_tukey_'+bin_labels[i], data=np.array(sat_nstar_tukey[i]))
            f.create_dataset('sat_large_scale_'+bin_labels[i], data=np.array(sat_nstar_large_scale[i]))
            f.create_dataset('sat_small_scale_'+bin_labels[i], data=np.array(sat_nstar_small_scale[i]))




