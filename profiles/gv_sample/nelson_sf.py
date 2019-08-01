import caesar
import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
import os
import matplotlib.colors as colors

plt.rc('text', usetex=True)
plt.rc('font', family='serif', )
plt.rcParams.update({'font.size': 14})

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
        new_cmap = colors.LinearSegmentedColormap.from_list(
                                'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
                                cmap(np.linspace(minval, maxval, n)))
        return new_cmap

def ssfr_line(masses, ssfr=1e-10):
    tens = np.array([10]*len(masses))
    return np.log10(np.power(tens, masses)*ssfr)

cmap = plt.get_cmap('jet_r')
new_cmap = truncate_colormap(cmap, 0.15, 0.95)

masses = np.arange(9.0, 13.0, 0.5)
line = ssfr_line(np.array(masses))

model = 'm100n1024'
wind = 's50j7k'
snap = '105'

h5_dir = '/home/sapple/simba_sizes/profiles/gv_sample/high_redshift/'
if not os.path.exists(h5_dir):
        os.makedirs(h5_dir)

plots_dir = h5_dir +model+'_'+snap+'/'
if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)

data_dir = '/home/rad/data/'+model+'/'+wind+'/'

sim =  caesar.load(data_dir+'Groups/'+model+'_'+snap+'.hdf5', LoadHalo=False)

gal_cent = np.array([i.central for i in sim.galaxies])
gal_sm = np.array([i.masses['stellar'].in_units('Msun') for i in sim.galaxies])
gal_sfr = np.array([i.sfr.in_units('Msun/yr') for i in sim.galaxies])
gal_ssfr = np.log10(gal_sfr / gal_sm)

gal_sm = np.log10(gal_sm)
gal_sfr = np.log10(gal_sfr)

sf_mask = gal_ssfr > 10.
with h5py.File(h5_dir+'sf_samples.h5', 'a') as f:
        f.create_dataset(model+'_'+snap, data=np.array(sf_mask))

plt.plot(masses, line, ls='--', lw=1.5, c='m', label=r'$\textrm{sSFR} = 10^{-10} \textrm{yr}^{-1}$')
plt.axvline(10., ls='--', lw=1.5, c='k')
plt.axvline(10.5, ls='--', lw=1.5, c='k')
plt.axvline(11., ls='--', lw=1.5, c='k')
plt.scatter(gal_sm, gal_sfr, s=1, c=gal_ssfr, cmap=new_cmap)
plt.xlim(9.5,12.5)
plt.ylim(-3.5, )
plt.clim(-12, -8)
plt.colorbar(label=r'$\textrm{log} (\textrm{sSFR} / \textrm{yr}^{-1})$')
plt.legend(fontsize=12)
plt.xlabel(r'$\log\ (M_{*} / M_{\odot})$')
plt.ylabel(r'$\textrm{log} (\textrm{SFR} / M_{\odot}\textrm{yr}^{-1})$')
plt.savefig(plots_dir+'sample_all_colormap.png')
plt.show()
plt.clf()


