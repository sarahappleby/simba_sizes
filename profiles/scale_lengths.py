import numpy as np
import h5py
import caesar
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

def exp(x, a, b):
    return a*np.exp(-1.*x/b)

def linear(x, a, b):
    return a*x + b

cent_sat = ['centrals', 'satellites']
ssfr_selection = ['star_forming', 'green_valley']
masks = [2, 3, 4]

save_dir = '/home/sapple/simba_sizes/profiles/paper/scale_lengths/'
plot_dir = '/home/sapple/simba_sizes/profiles/paper/scale_lengths/plots/'

model = 'm100n1024'
wind = 's50j7k'
snap = '151'
start_r = 1.

data_dir = '/home/rad/data/'+model+'/'+wind+'/Groups/'
sim = caesar.load(data_dir+model+'_'+snap+'.hdf5', LoadHalo=False)

gal_sm = np.array([i.masses['stellar'].in_units('Msun') for i in sim.galaxies])
gal_sfr = np.array([i.sfr.in_units('Msun/yr') for i in sim.galaxies])
gal_ssfr = gal_sfr / gal_sm

halflight_file = '/home/sapple/simba_sizes/sizes/data/halfradius_R.h5'
with h5py.File(halflight_file, 'r') as f:
    gal_rad = f[model+'_'+wind+'_'+snap+'_halflight'][:] # these are in pkpc
gal_rad = np.sum(gal_rad, axis=0) / 3.

h1_scalelengths = np.zeros(len(sim.galaxies))
sfr_scalelengths = np.zeros(len(sim.galaxies))

for select in cent_sat:
    
    print('starting ' + select)

    for ssfr_s in ssfr_selection:

        print('starting '+ ssfr_s)

        profile_dir = '/home/sapple/simba_sizes/profiles/paper/'+select+'/m100n1024_151/s50j7k/'+ssfr_s+'/rotated_faceon/'
        for m in masks:
            with h5py.File(profile_dir+'/mask_'+str(m)+'_all_profiles.h5', 'r') as f:
                gal_ids = f['gal_ids'].value
                h1_profs = f['h1_m'].value
                sfr_profs = f['gas_sfr'].value

            dr = 0.2
            n = h1_profs.shape[1]
            factor = dr*n
            bins = np.arange(0., factor, dr)
            rplot = bins+(dr*0.5)
            start = int(start_r / dr)

            for i, prof in enumerate(h1_profs):
                rad = gal_rad[gal_ids[i]]
                prof[np.where(prof == 0.)[0]] = 1.
                try:
                    popt, pcov = curve_fit(linear, (rad*rplot)[start:], np.log(prof)[start:])
                    plt.plot(rplot*rad, np.log(prof))
                    y = linear((rplot*rad)[start:], popt[0], popt[1])
                    plt.plot((rplot*rad)[start:], y)
                    plt.savefig(plot_dir+'h1_gal_'+str(gal_ids[i])+'.png')
                    plt.clf()
                    
                    h1_scalelengths[gal_ids[i]] = -1.0 / popt[0]


                except RuntimeError:
                    h1_scalelengths[gal_ids[i]] = 0.

            for i, prof in enumerate(sfr_profs):
                rad = gal_rad[gal_ids[i]]
                prof[np.where(prof == 0.)[0]] = 1.e-6
                try:
                    popt, pcov = curve_fit(linear, (rad*rplot)[start:], np.log(prof)[start:])
                    plt.plot(rplot*rad, np.log(prof))
                    y = linear((rplot*rad)[start:], popt[0], popt[1])
                    plt.plot((rplot*rad)[start:], y)
                    plt.savefig(plot_dir+'sfr_gal_'+str(gal_ids[i])+'.png')
                    plt.clf()

                    sfr_scalelengths[gal_ids[i]] = -1.0 / popt[0]

                except RuntimeError:
                    sfr_scalelengths[gal_ids[i]] = 0.

with h5py.File(save_dir+'scale_lengths.h5', 'a') as f:
    f.create_dataset('h1_scale_length', data=np.array(h1_scalelengths))
    f.create_dataset('sfr_scale_length', data=np.array(sfr_scalelengths))
    f.create_dataset('sm', data=np.array(gal_sm))
    f.create_dataset('ssfr', data=np.array(gal_ssfr))


