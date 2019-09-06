import numpy as np
import h5py
import caesar

def broeils_rhee(log_mh1):
    return (0.51 *log_mh1) - 3.33 

model = 'm100n1024'
wind = 's50j7k'
snap = '151'
mass_bins = [10., 10.5, 11.]

if wind == 's50j7k':
    halflight_file = '/home/sapple/simba_sizes/sizes/data/halfradius_R.h5'
else:
    halflight_file = '/home/sapple/simba_sizes/sizes/data/halfradius_agn_R.h5'

with h5py.File(halflight_file, 'r') as f:
    gal_rad = f[model+'_'+wind+'_'+snap+'_halflight'][:] # these are in pkpc
gal_rad = np.sum(gal_rad, axis=0) / 3.


data_dir = '/home/rad/data/'+model+'/'+wind+'/'
sim =  caesar.load(data_dir+'Groups/'+model+'_'+snap+'.hdf5', LoadHalo=False)
h = sim.simulation.hubble_constant
redshift = sim.simulation.redshift

gal_sm = np.array([i.masses['stellar'].in_units('Msun') for i in sim.galaxies])
gal_sfr = np.array([i.sfr.in_units('Msun/yr') for i in sim.galaxies])
gal_ssfr = gal_sfr / gal_sm
gal_h1 = np.array([i.masses['HI'].in_units('Msun') for i in sim.galaxies])
gal_h2 = np.array([i.masses['H2'].in_units('Msun') for i in sim.galaxies])

gal_sm = np.log10(gal_sm)
gal_ssfr = np.log10(gal_ssfr)
gal_h1 = np.log10(gal_h1)
gal_h2 = np.log10(gal_h2)

sf_sample_file = '/home/sapple/simba_sizes/profiles/gv_sample/belfiore/s50j7k/selection_2/sf_samples.h5'
with h5py.File(sf_sample_file, 'r') as f:
    sf_gals = f[model+'_'+snap].value

gv_sample_file = '/home/sapple/simba_sizes/profiles/gv_sample/belfiore/s50j7k/selection_2/gv_samples.h5'
with h5py.File(gv_sample_file, 'r') as f:
    gv_gals = f[model+'_'+snap].value

print('Data for star forming galaxies: ') 
for j in range(len(mass_bins)):
        print('\n')
        print('Looking at mass bin ' +str(j) )
        if j != len(mass_bins) - 1:
                sm_mask = (gal_sm > mass_bins[j]) & (gal_sm < mass_bins[j+1])
        else:
                sm_mask = gal_sm > mass_bins[j]

        mask = sm_mask * sf_gals
        print('Median M* : ' + str(np.median(gal_sm[mask])))
        print('Median SFR : ' + str(np.median(gal_sfr[mask])))
        print('Median sSFR : ' + str(np.median(gal_ssfr[mask])))
        print('Median HI : ' + str(np.median(gal_h1[mask])))
        print('Median H2 : ' + str(np.median(gal_h2[mask])))

        med_rhalf = np.median(gal_rad[mask])
        med_dhi = 10.**(broeils_rhee(np.median(gal_h1[mask])))
        print('Median Rhalf : ' + str(med_rhalf))
        print('Median DHI : ' + str(med_dhi))
        print('Ratio : ' +str(med_dhi / med_rhalf))
        print('\n')


print('Data for green valley galaxies: ')
for j in range(len(mass_bins)):
        print('\n')
        print('Looking at mass bin ' +str(j) )
        if j != len(mass_bins) - 1:
                sm_mask = (gal_sm > mass_bins[j]) & (gal_sm < mass_bins[j+1])
        else:
                sm_mask = gal_sm > mass_bins[j]

        mask = sm_mask * gv_gals
        print('Median M* : ' + str(np.median(gal_sm[mask])))
        print('Median SFR : ' + str(np.median(gal_sfr[mask])))
        print('Median sSFR : ' + str(np.median(gal_ssfr[mask])))
        print('Median HI : ' + str(np.median(gal_h1[mask])))
        print('Median H2 : ' + str(np.median(gal_h2[mask])))
        
        med_rhalf = np.median(gal_rad[mask])
        med_dhi = 10.**(broeils_rhee(np.median(gal_h1[mask])))

        print('Median Rhalf : ' + str(med_rhalf))
        print('Median DHI : ' + str(med_dhi))
        print('Ratio : ' +str(med_dhi / med_rhalf))


        print('\n')

