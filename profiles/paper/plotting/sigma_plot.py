# plot surface densities of stellar mass against star formation rate for all profiles?? 
# start off just for star forming galaxies, we can see where the krumholz line lies here

# read in data, get the half light radii, scale each profile by radius to get colour scale

import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import caesar
import sys
from plotting_methods import plot_kennicutt

sys.path.append('/home/sapple/tools')
import plotmedian as pm

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
        new_cmap = colors.LinearSegmentedColormap.from_list(
                                'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
                                cmap(np.linspace(minval, maxval, n)))
        return new_cmap

selection = sys.argv[1]

model = 'm100n1024'
wind = 's50j7k'
snap = '151'
sample_file = '/home/sapple/simba_sizes/profiles/gv_sample/belfiore/s50j7k/selection_2/'+selection+'_samples.h5'
angle = 'rotated_faceon'
mlim = 10.

cmap = plt.get_cmap('jet')
new_cmap = truncate_colormap(cmap, 0.1, 1.)


if selection == 'sf':
    name = 'star_forming'
elif selection == 'gv':
    name = 'green_valley'

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 14})


with h5py.File(sample_file, 'r') as f:
    id_mask = f[model+'_'+snap].value

halflight_file = '/home/sapple/simba_sizes/sizes/data/halfradius_R.h5'
with h5py.File(halflight_file, 'r') as f:
    gal_rad = f[model+'_'+wind+'_'+snap+'_halflight'][:] # these are in pkpc
gal_rad = np.sum(gal_rad, axis=0) / 3.

data_dir = '/home/rad/data/'+model+'/'+wind+'/'
sim =  caesar.load(data_dir+'Groups/'+model+'_'+snap+'.hdf5', LoadHalo=False)
gal_sm = np.log10(np.array([i.masses['stellar'].in_units('Msun') for i in sim.galaxies]))
mask = (gal_sm > mlim) & id_mask

basic_dir = '/home/sapple/simba_sizes/profiles/paper/'
plot_dir = basic_dir+'plotting/plots/'
profs_file = basic_dir + 'all_profs/'+name+'_'+angle + '.h5'

gal_rad = gal_rad[mask]

with h5py.File(profs_file, 'r') as f:
    h1_h2_m = np.log10(f['h2_m'].value[mask] + f['h1_m'].value[mask]) - 6. # go from kpc to pc
    sfr = np.log10(f['sfr'].value[mask] + 1.e-6)

n = sfr.shape[1]
dr = 0.2
factor = dr*n
bins = np.arange(0., factor, dr)
rplot = bins+(dr*0.5)
rplot = np.tile(rplot, len(sfr))

h1_h2_m = h1_h2_m.flatten()
sfr = sfr.flatten()
bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(h1_h2_m[h1_h2_m > -1], sfr[h1_h2_m > -1], stat='mean')

fig, ax = plt.subplots()
plot_kennicutt(ax, labels=['Kennicutt 1998', 'Bigiel 2008'])
ax.plot(bin_cent, ymean, c='m', lw=2.5, ls='--')
im = ax.scatter(h1_h2_m.flatten(), sfr.flatten(), c=rplot, s=0.5, cmap =new_cmap)
cbar = fig.colorbar(im,ax=ax, label=r'$ R / R_{half}$')
ax.set_xlabel(r'$ \textrm{log} (\Sigma _{HI + H2}/ M_{\odot}\textrm{pc}^{-2}) $')
ax.set_ylabel(r'$ \textrm{log} (\Sigma _{\textrm{SFR}}/ M_{\odot}\textrm{yr}^{-1} \textrm{kpc}^{-2}) $')
cbar.set_clim(0, 5)
ax.set_xlim(-1, 4)
ax.set_ylim(-7,)
ax.legend()
plt.savefig(plot_dir+'sigma_sfr_gas_'+selection+'.png')
plt.clf()
