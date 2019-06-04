import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np 
import caesar
import h5py
import os
import sys
sys.path.append('/home/sapple/tools')
import plotmedian as pm

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
    new_cmap = colors.LinearSegmentedColormap.from_list('trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
                                                        cmap(np.linspace(minval, maxval, n)))
    return new_cmap

def concentration(r20, r80):
    return 5.*np.log10(r80/r20)

radius_file = '/home/sapple/simba_sizes/concentration/data/r20r80.h5'
plot_dir = '/home/sapple/simba_sizes/concentration/plots/'

model = 'm100n1024'
wind = 's50'
snaps = ['062', '078', '090', '105', '125', '151']
simba_z = [3.0, 2.0, 1.5, 1.0, 0.5, 0.0]

mmin = 8.7
mmax = 12.5

cmap = plt.get_cmap('jet_r')
cmap = truncate_colormap(cmap, 0.05, 1.0)

fig = plt.figure(figsize=(15, 17))
for i, snap in enumerate(snaps):

    caesar_data = '/home/sapple/simba_sizes/sizes/data/'+model+'_'+wind+'_'+snap+'_caesar_data.h5'
    if not os.path.isfile(caesar_data):
        infile = '/home/rad/data/'+model+'/'+wind+'/Groups/'+model+'_'+snap+'.hdf5'
        sim = caesar.load(infile, LoadHalo=False)
        central = np.asarray([j.central for j in sim.galaxies])
        ms = np.asarray([j.masses['stellar'].in_units('Msun') for j in sim.galaxies])
        sfr = np.array([j.sfr.in_units('Msun/yr') for j in sim.galaxies])
        sfr[np.where(sfr == 1.)[0]] = 0.
        with h5py.File(caesar_data, 'w') as f:
            f.create_dataset('central', data=np.array(central))
            f.create_dataset('stellar_mass', data=np.array(ms))
            f.create_dataset('sfr', data=np.array(sfr))
    else:
        with h5py.File(caesar_data, 'r') as f:
            central = f['central'].value
            ms = f['stellar_mass'].value
            sfr = f['sfr'].value

    ssfr = np.log10(1.e9*sfr/ms)
    ssfrlim = -1.8+0.3*simba_z[i]

    with h5py.File(radius_file, 'r') as r:
        r20 = r[model+'_'+wind+'_'+snap+'_r20'].value
        r80 = r[model+'_'+wind+'_'+snap+'_r80'].value
    r20 = np.sum(r20, axis=0) / 3.
    r80 = np.sum(r80, axis=0) / 3.
    conc = concentration(r20, r80)

    simba_c = ssfr
    simba_c[simba_c < -2.5] = -2.5
    pixsize = 1*(np.log10(ms)-min(np.log10(ms)))+0.5

    ax = fig.add_subplot(3, 2, i+1)
    im = ax.scatter(np.log10(ms), conc, c=simba_c, s=pixsize, lw=0, cmap=cmap, label='Simba')
    cbar = fig.colorbar(im,ax=ax, label=r'$\textrm{log} (\textrm{sSFR} / \textrm{Gyr}^{-1})$')
    cbar.ax.tick_params(labelsize=10)

    if simba_z[i] <= 1.5:
        bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(np.log10(ms[ssfr>ssfrlim]),conc[ssfr>ssfrlim])
        ax.plot(bin_cent,ymean,'--',lw=4,color='c', label='Star forming')
        bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(np.log10(ms[ssfr<ssfrlim]),conc[ssfr<ssfrlim])
        ax.plot(bin_cent,ymean,'--',lw=4,color='m', label='Passive')

    else:
        bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(np.log10(ms),conc)
        ax.plot(bin_cent,ymean,'--',lw=4,color='b', label='Simba')

    ax.annotate('z=%g'%(np.round(simba_z[i],1)), xy=(0.1, 0.9), xycoords='axes fraction',size=16,bbox=dict(boxstyle="round", fc="w"))
    ax.minorticks_on()
    ax.set_xlim(mmin,mmax)
    ax.set_xlabel(r'$\log\ (M_{*} / M_{\odot})$',fontsize=16)
    ax.set_ylabel(r'$C$')
    ax.legend(loc='lower right', fontsize=8)

plt.savefig(plot_dir+'conc_'+model+'.png', bbox_inches='tight', format='png')
plt.clf()

