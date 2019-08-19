import numpy as np
import h5py
import caesar
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def exp(x, a, b):
    return a*np.exp(-1.*x/b)

cent_sat = ['centrals', 'satellites']
ssfr_selection = ['star_forming', 'green_valley']
masks = [1, 2, 3, 4]

save_dir = '/home/sapple/simba_sizes/profiles/paper/scale_lengths/'
plot_dir = '/home/sapple/simba_sizes/profiles/paper/scale_lengths/plots/'

model = 'm100n1024'
wind = 's50j7k'
snap = '151'
start_r = 0.8

data_dir = '/home/rad/data/'+model+'/'+wind+'/Groups/'
sim = caesar.load(data_dir+model+'_'+snap+'.hdf5', LoadHalo=False)

gal_sm = np.array([i.masses['stellar'].in_units('Msun') for i in sim.galaxies])
gal_sfr = np.array([i.sfr.in_units('Msun/yr') for i in sim.galaxies])
gal_ssfr = gal_sfr / gal_sm

halflight_file = '/home/sapple/simba_sizes/sizes/data/halfradius_R.h5'
with h5py.File(halflight_file, 'r') as f:
    gal_rad = f[model+'_'+wind+'_'+snap+'_halflight'][:] # these are in pkpc
gal_rad = np.sum(gal_rad, axis=0) / 3.

scale_lengths = np.zeros(len(sim.galaxies))

for select in cent_sat:
    
    for ssfr_s in ssfr_selection:

        profile_dir = '/home/sapple/simba_sizes/profiles/paper/'+select+'/m100n1024_151/s50j7k/'+ssfr_s+'/random_orientation/'
        for m in masks:
            with h5py.File(profile_dir+'/mask_'+str(m)+'_all_profiles.h5', 'r') as f:
                gal_ids = f['gal_ids'].value
                h1_profs = f['h1'].value
        
            dr = 0.2
            n = h1_profs.shape[1]
            factor = dr*n
            bins = np.arange(0., factor, dr)
            rplot = bins+(dr*0.5)
            start = int(start_r / dr)

            for i, prof in enumerate(h1_profs):
                rad = gal_rad[gal_ids[i]]
                popt, pcov = curve_fit(exp, (rad*rplot)[start:], prof[start:])
                scale_lengths[gal_ids[i]] = popt[1]

                plt.plot(rplot*rad, prof)
                x = exp((rplot*rad)[:start], popt[0], popt[1])
                plt.plot((rplot*rad)[:start], x)
                plt.savefig(plot_dir+'gal_'+str(gal_ids[i])+'.png')
                plt.clf()

with h5py.File(save_dir+'scale_lengths.h5', 'a') as f:
    f.create_dataset('scale_length', data=np.array(scale_lengths))
    f.create_dataset('sm', data=np.array(gal_sm))
    f.create_dataset('ssfr', data=np.array(gal_ssfr))
