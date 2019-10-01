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

caesar_dir = '/home/rad/data/'+model+'/'+wind+'/Groups/'
sim =  caesar.load(caesar_dir+model+'_'+snap+'.hdf5', LoadHalo=False)
gal_pos = np.array([i.pos.in_units('kpc/h') for i in sim.galaxies])
boxsize = sim.simulation.boxsize.in_units('kpc/h')

for i, m in enumerate(masks):

    with h5py.File(centrals_dir+'mask_'+str(m)+'_all_profiles.h5', 'r') as f:
        cen_gas_h2_z = f['h2-weighted_z'].value
        cen_gal_ids = f['gal_ids'].value

    cen_pos = gal_pos[cen_gal_ids]

    if gals == 'all':
        with h5py.File(sats_dir+'mask_'+str(m)+'_all_profiles.h5', 'r') as f:
            sat_gas_h2_z = f['h2-weighted_z'].value
            sat_gal_ids = f['gal_ids'].value

        sat_pos = gal_pos[sat_gal_ids]

    if i == 0:
        n = cen_gas_h2_z.shape[1]
        cen_h2_z_jk = np.zeros((len(masks), n)); cen_h2_z_cv_err = np.zeros((len(masks), n))

        if gals == 'all':

            sat_h2_z_jk = np.zeros((len(masks), n)); sat_h2_z_cv_err = np.zeros((len(masks), n))

            all_h2_z_jk = np.zeros((len(masks), n)); all_h2_z_cv_err = np.zeros((len(masks), n))

    cen_h2_z_jk[i], cen_h2_z_cv_err[i] = cosmic_variance(cen_gas_h2_z, cen_pos, boxsize, 'h2_z')

    if gals == 'all':

        sat_h2_z_jk[i], sat_h2_z_cv_err[i] = cosmic_variance(sat_gas_h2_z, sat_pos, boxsize, 'h2_z')

        h2_z = np.concatenate((cen_gas_h2_z, sat_gas_h2_z))
        all_pos = np.concatenate((cen_pos, sat_pos))

        all_h2_z_jk[i], all_h2_z_cv_err[i] = cosmic_variance(h2_z, all_pos, boxsize, 'h2_z')

    with h5py.File(results_dir+'_h2-weighted_z_data.h5', 'a') as f: 
        f.create_dataset('cen_jk_'+bin_labels[i], data=np.array(cen_h2_z_jk[i]))
        f.create_dataset('cen_cv_err_'+bin_labels[i], data=np.array(cen_h2_z_cv_err[i]))
        if gals == 'all':
            f.create_dataset('sat_jk_'+bin_labels[i], data=np.array(sat_h2_z_jk[i]))
            f.create_dataset('sat_cv_err_'+bin_labels[i], data=np.array(sat_h2_z_cv_err[i]))
            f.create_dataset('all_jk_'+bin_labels[i], data=np.array(all_h2_z_jk[i]))
            f.create_dataset('all_cv_err_'+bin_labels[i], data=np.array(all_h2_z_cv_err[i]))


