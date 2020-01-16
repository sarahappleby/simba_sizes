# atm i have the profiles saved in different mass bins, labelled masks 2, 3, 4 for silly reasons
# rather than running again (stupid) I should just join these bins and take slightly different mass cuts

# read in the data from each mask, get the gal ids, assign profiles to the place in an empty array
# make for centrals, satellites, both

import h5py
import numpy as np
import sys
import caesar

model = sys.argv[1]
wind = sys.argv[2]
snap = sys.argv[3]
selection = sys.argv[4]
angle = sys.argv[5]

if selection =='gv':
    name = 'green_valley'
elif selection == 'sf':
    name = 'star_forming'

basic_dir = '/home/sapple/simba_sizes/profiles/paper/'
centrals_dir = basic_dir + 'centrals/'+model+'_'+snap+'/'+wind+'/'+name+'/'+angle+'/'
sats_dir = basic_dir + 'satellites/'+model+'_'+snap+'/'+wind+'/'+name+'/'+angle+'/'
results_dir = basic_dir+'all_profs/'

data_dir = '/home/rad/data/'+model+'/'+wind+'/'
sim =  caesar.load(data_dir+'Groups/'+model+'_'+snap+'.hdf5', LoadHalo=False)

masks = [2, 3, 4]
n = 25

all_star_m = np.zeros((len(sim.galaxies), n))
all_h1 = np.zeros((len(sim.galaxies), n))
all_h2 = np.zeros((len(sim.galaxies), n))
all_h1_m = np.zeros((len(sim.galaxies), n))
all_h2_m = np.zeros((len(sim.galaxies), n))
all_gm = np.zeros((len(sim.galaxies), n))
all_sfr = np.zeros((len(sim.galaxies), n))
all_ages = np.zeros((len(sim.galaxies), n))
all_temp = np.zeros((len(sim.galaxies), n))

for i, m in enumerate(masks):
        with h5py.File(centrals_dir+'mask_'+str(m)+'_all_profiles.h5', 'r') as f:
            gal_ids = f['gal_ids'].value
            all_star_m[gal_ids] = f['sm'].value
            all_sfr[gal_ids] = f['gas_sfr'].value
            all_h1[gal_ids] = f['h1'].value
            all_h2[gal_ids] = f['h2'].value
            all_h1_m[gal_ids] = f['h1_m'].value
            all_h2_m[gal_ids] = f['h2_m'].value
            all_temp[gal_ids] = f['temp'].value
            all_ages[gal_ids] = f['ages'].value
            all_gm[gal_ids] = f['gm'].value

        with h5py.File(sats_dir+'mask_'+str(m)+'_all_profiles.h5', 'r') as f:
            gal_ids = f['gal_ids'].value
            all_star_m[gal_ids] = f['sm'].value
            all_sfr[gal_ids] = f['gas_sfr'].value
            all_h1[gal_ids] = f['h1'].value
            all_h2[gal_ids] = f['h2'].value
            all_h1_m[gal_ids] = f['h1_m'].value
            all_h2_m[gal_ids] = f['h2_m'].value
            all_temp[gal_ids] = f['temp'].value
            all_ages[gal_ids] = f['ages'].value
            all_sfr[gal_ids] = f['gas_sfr'].value

with h5py.File(results_dir+name+'_'+angle+'.h5', 'a') as f:
    f.create_dataset('sm', data=all_star_m)
    f.create_dataset('h1', data=all_h1)
    f.create_dataset('h2', data=all_h2)
    f.create_dataset('h1_m', data=all_h1_m)
    f.create_dataset('h2_m', data=all_h2_m)
    f.create_dataset('gm', data=all_gm)
    f.create_dataset('sfr', data=all_sfr)
    f.create_dataset('temp', data=all_temp)
    f.create_dataset('ages', data=all_ages)
