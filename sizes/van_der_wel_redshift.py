import matplotlib
import matplotlib.pyplot as plt
import os
import numpy as np
import h5py
from scipy.optimize import curve_fit

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def evolution(z, a, b):
    return a*((1+z)**b)

radius_file = '/home/sapple/simba_sizes/sizes/data/halfradius_V.h5'
plot_dir = '/home/sapple/simba_sizes/sizes/plots/'
simba_caesar = '/home/sapple/simba_sizes/sizes/data/'

model = 'm100n1024'
wind = 's50'

m_star = 5.e10
m_star_str = str(m_star)[0] + 'e' + str(int(np.log10(m_star)))

simba_snaps = ['062', '078', '090', '105', '125', '145', '151']
simba_z = [3.0, 2.0, 1.5, 1.0, 0.5, 0.1, 0.0]

vdw_blue = [8.9 * (1+i)**-0.75 for i in np.arange(0., 3.1, 0.1)]
vdw_red = [5.6* (1+i)**-1.48 for i in np.arange(0., 3.1, 0.1)]

blue_sizes = np.zeros(len(simba_z))
blue_per_25 = np.zeros(len(simba_z))
blue_per_75 = np.zeros(len(simba_z))

red_sizes = np.zeros(len(simba_z))
red_per_25 = np.zeros(len(simba_z))
red_per_75 = np.zeros(len(simba_z))

for i, snap in enumerate(simba_snaps):
    with h5py.File(radius_file, 'r') as r:
        rhalf = r[model+'_'+wind+'_'+snap+'_halflight'][:]
    rhalf = np.sum(rhalf, axis=0) / 3.
    caesar_data = simba_caesar + model+'_'+wind+'_'+snap+'_caesar_data.h5'
    if not os.path.isfile(caesar_data):
        import caesar
        infile = '/home/rad/data/'+model+'/'+wind+'/Groups/'+model+'_'+snap+'.hdf5'
        sim = caesar.load(infile, LoadHalo=False)
        central = np.asarray([j.central for j in sim.galaxies])
        smass = np.asarray([j.masses['stellar'].in_units('Msun') for j in sim.galaxies])
        sfr = np.array([j.sfr.in_units('Msun/yr') for j in sim.galaxies])
        sfr[np.where(sfr == 1.)[0]] = 0.
        with h5py.File(caesar_data, 'w') as f:
            f.create_dataset('central', data=np.array(central))
            f.create_dataset('stellar_mass', data=np.array(smass))
            f.create_dataset('sfr', data=np.array(sfr))
    else:
        with h5py.File(caesar_data, 'r') as f:
            central = f['central'][:]
            smass = f['stellar_mass'][:]
            sfr = f['sfr'][:]
    ssfr = np.log10(1.e9*sfr/smass)
    ssfrlim = -1.8+0.3*simba_z[i]

    mask = (smass > 0.9*m_star) & (smass < 1.1*m_star)
    star_forming = ssfr > ssfrlim

    blue_sizes[i] = np.median(rhalf[central*mask*star_forming])
    blue_per_25[i] = np.percentile(rhalf[central*mask*star_forming], 25)
    blue_per_75[i] = np.percentile(rhalf[central*mask*star_forming], 75)

    red_sizes[i] = np.median(rhalf[central*mask*np.invert(star_forming)])
    red_per_25[i] = np.percentile(rhalf[central*mask*np.invert(star_forming)], 25)
    red_per_75[i] = np.percentile(rhalf[central*mask*np.invert(star_forming)], 75)

    print 'z = ' + str(simba_z[i])
    print 'Number of star forming galaxies: ' + str(len(rhalf[central*mask*star_forming]))
    print 'Number of passive galaxies: ' + str(len(rhalf[central*mask*np.invert(star_forming)]))
    print '\n'

popt, pcov = curve_fit(evolution, simba_z, blue_sizes)
blue_string = r'$R/\textrm{kpc} = %.1f (1 + z)^{%.2f}$' % (popt[0], popt[1])
blue_fit = evolution(np.array(simba_z), popt[0], popt[1])

popt, pcov = curve_fit(evolution, simba_z, red_sizes)
red_string = r'$R/\textrm{kpc} = %.1f (1 + z)^{%.2f}$' % (popt[0], popt[1])
red_fit = evolution(np.array(simba_z), popt[0], popt[1])

plt.plot(0., 0.81, color='b', marker='s', markersize=5., linestyle='None', label='SDSS-SF')
plt.plot(0., 0.69, color='r', marker='s', markersize=5., linestyle='None', label='SDSS-Q')

# plotting the medians with the data
plt.plot(np.arange(0., 3.1, 0.1), np.log10(vdw_blue), '-', lw=1.5, color='b', label='CANDELS-LTG (van der Wel+14)')
plt.plot(np.arange(0., 3.1, 0.1), np.log10(vdw_red), '-', lw=1.5, color='r', label='CANDELS-ETG (van der Wel+14)')

plt.plot(simba_z, np.log10(blue_sizes), linestyle='--', c='c', lw=1.5, label='Star forming; '+ blue_string)
plt.fill_between(simba_z, np.log10(blue_per_25), np.log10(blue_per_75), facecolor='c', alpha=0.15, linewidth=1)

plt.plot(simba_z, np.log10(red_sizes), linestyle='--', c='m', lw=1.5, label='Passive; '+red_string)
plt.fill_between(simba_z, np.log10(red_per_25), np.log10(red_per_75), facecolor='m', alpha=0.15, linewidth=1)

plt.xlabel(r'z', fontsize=16)
#plt.xlim(3.1, -0.1)
plt.ylabel(r'$\textrm{log}\ ( R_\textrm{half}\ /\ \textrm{kpc})$' ,fontsize=16)
plt.legend(loc=1)
plt.savefig(plot_dir+'redshift_medians_'+m_star_str+'.png')
plt.show()
plt.clf()
