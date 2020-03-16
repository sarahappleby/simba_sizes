import numpy as np
import h5py
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import caesar
import sys
from plotting_methods import plot_kennicutt
from bigiel_data import *

sys.path.append('/home/sapple/tools')
import plotmedian as pm

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
        new_cmap = colors.LinearSegmentedColormap.from_list(
                                'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
                                cmap(np.linspace(minval, maxval, n)))
        return new_cmap

def get_area(dr, n):
    area = np.zeros(n)
    for i in range(n):
        if i == 0.:
            area[i] = np.pi * (dr**2.)
        else:
            area[i] = np.pi * ((dr*(i+1))**2 - (dr*i)**2)
    return area


model = 'm100n1024'
wind = 's50'
snap = '151'
angle = 'random_orientation'
mlim = 10.

bigiel_to_salpeter = 1.59
salpeter_to_chabrier = 1.75
bigiel_imf_factor = bigiel_to_salpeter / salpeter_to_chabrier

m51_kennicutt_07 = np.transpose(np.array(m51_kennicutt_07))
spirals_kennicutt_98 = np.transpose(np.array(spirals_kennicutt_98))
bigiel_contour_1 = np.transpose(np.array(bigiel_contour_1))
bigiel_contour_2 = np.transpose(np.array(bigiel_contour_2))
bigiel_contour_3 = np.transpose(np.array(bigiel_contour_3))
bigiel_contour_4 = np.transpose(np.array(bigiel_contour_4))

cmap = plt.get_cmap('jet_r')
new_cmap = truncate_colormap(cmap, 0.1, 1.)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 12})

halflight_file = '/home/sapple/simba_sizes/sizes/data/pyloser_sdss_r.h5'
with h5py.File(halflight_file, 'r') as f:
    gal_rad = f['abs_'+model+'_'+wind+'_'+snap][:] # these are in pkpc
gal_rad = np.sum(gal_rad, axis=0) / 3.

data_dir = '/home/rad/data/'+model+'/'+wind+'/'
sim =  caesar.load(data_dir+'Groups/'+model+'_'+snap+'.hdf5', LoadHalo=False)
redshift = sim.simulation.redshift

gal_sm = np.array([i.masses['stellar'].in_units('Msun') for i in sim.galaxies])
gal_sfr = np.array([i.sfr.in_units('Msun/Gyr') for i in sim.galaxies])
gal_ssfr = np.log10((gal_sfr / gal_sm) + 10**(-2.5))
gal_sm = np.log10(gal_sm)

mass_mask = (gal_sm > mlim)
ssfr_mask = gal_ssfr > (-1.8 + 0.3 * redshift)

choose_mask = mass_mask * ssfr_mask
gal_rad = gal_rad[choose_mask]
gal_ssfr = gal_ssfr[choose_mask]

basic_dir = '/home/sapple/simba_sizes/profiles/paper/'
plot_dir = basic_dir+'plotting/plots/'
profs_file = basic_dir + '/with_dust_sizes/'+wind+'/'+model+'_'+snap+'/all_profiles_'+angle+'.h5'

with h5py.File(profs_file, 'r') as f:
    h1_surf_dens = f['h1_m'][:][choose_mask]
    h2_surf_dens = f['h2_m'][:][choose_mask]
    sfr_surf_dens = f['gas_sfr'][:][choose_mask]

n = 5
dr = 0.2

gal_rad_tile = np.transpose(np.tile(gal_rad, (n, 1)))
area = get_area(dr, n)
area_phys = gal_rad_tile**2.  * area 
area_total = np.sum(area_phys[:, :n], axis=1)
h1_mass = h1_surf_dens[:, :n] * area_phys
h2_mass = h2_surf_dens[:, :n] * area_phys
sfr = sfr_surf_dens[:, :n] * area_phys

gas_mass_total = np.sum(h1_mass, axis=1) + np.sum(h2_mass, axis=1)
sigma_h2 = np.log10(np.sum(h2_mass, axis=1) / area_total) - 6.
sigma_gas = np.log10(gas_mass_total / area_total) - 6.
sfr_total = np.sum(sfr, axis=1)
sigma_sfr = np.log10((sfr_total / area_total) + 10**-4.5) 

condition = (sigma_gas > 0.)

bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(sigma_gas[condition], sigma_sfr[condition], stat='mean')

fig, ax = plt.subplots(figsize=(5, 4))

if snap == '151':
    plot_kennicutt(ax, labels=[r'$\Sigma_{\textrm{SFR}} \propto \Sigma_{\textrm{HI} + \textrm{H}_2}^{1.4}$'])
    ax.plot(bigiel_contour_1[0], bigiel_contour_1[1] * bigiel_imf_factor, c='k', lw=1, label='Bigiel+08')
    ax.plot(bigiel_contour_2[0], bigiel_contour_2[1] * bigiel_imf_factor, c='k', lw=1)
    ax.plot(bigiel_contour_3[0], bigiel_contour_3[1] * bigiel_imf_factor, c='k', lw=1)
    ax.plot(bigiel_contour_4[0], bigiel_contour_4[1] * bigiel_imf_factor, c='k', lw=1)
#ax.scatter(spirals_kennicutt_98[0], spirals_kennicutt_98[1] * bigiel_imf_factor, c='k', marker='o', s=10, label='Spirals, Kennicutt 98')
#ax.scatter(m51_kennicutt_07[0], m51_kennicutt_07[1], c='k', marker='+', s=14, label='M51, Kennicutt+07')

ax.plot(bin_cent, ymean, c='dodgerblue', lw=2.5, ls='--', label=r'$\textrm{HI} + \textrm{H}_2$')
im = ax.scatter(sigma_gas[gal_ssfr>-1.8], sigma_sfr[gal_ssfr> -1.8], c=gal_ssfr[gal_ssfr>-1.8], s=0.5, cmap =new_cmap)
cbar = fig.colorbar(im,ax=ax, label=r'$ \textrm{log} (sSFR / Gyr^{-1})$')
ax.set_xlabel(r'$ \textrm{log} (\Sigma _{HI + H2}/ M_{\odot}\textrm{pc}^{-2}) $')
ax.set_ylabel(r'$ \textrm{log} (\Sigma _{\textrm{SFR}}/ M_{\odot}\textrm{yr}^{-1} \textrm{kpc}^{-2}) $')

bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(sigma_h2[condition], sigma_sfr[condition], stat='mean')
ax.plot(bin_cent, ymean, c='m', lw=2.5, ls='--', label=r'$\textrm{H}_2$ only')

im.set_clim(-2, 1.)

if snap == '151':
    ax.set_xlim(-0.5, 5)
    ax.set_ylim(-5.,0.5)
else:
    ax.set_xlim(-0.5,5.)
    ax.set_ylim(-5.,)

ax.annotate('z=%g'%(np.round(sim.simulation.redshift,1)), xy=(0.05, 0.9), xycoords='axes fraction',size=16,bbox=dict(boxstyle="round", fc="w"))

ax.legend(loc=4, fontsize=10)
plt.savefig(plot_dir+'surface_densities_'+model+'_'+snap+'_'+angle+'.png')
