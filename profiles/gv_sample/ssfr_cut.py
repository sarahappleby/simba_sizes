import caesar
import numpy as np
import matplotlib.pyplot as plt
import h5py
import sys
import os
import matplotlib.colors as colors

def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=100):
        new_cmap = colors.LinearSegmentedColormap.from_list(
                                'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
                                cmap(np.linspace(minval, maxval, n)))
        return new_cmap

def sfms_line(x0, a=1., b=-10.):
	return x0*a + b

masses = [10., 10.5, 11.]
ssfr_min = -11.5
ssfr_max = -10.5
sm_plot = np.arange(9.0, 12., 0.5)

model = sys.argv[1]
wind = sys.argv[2]
snap = sys.argv[3]

h5_dir = '/home/sapple/simba_sizes/profiles/gv_sample/ssfr_cuts-10.5/'+wind+'/'
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
gal_sfr[np.where(gal_sfr == 1.)[0]] = 0.
gal_ssfr = np.log10(gal_sfr / gal_sm)

gal_sm = np.log10(gal_sm)
gal_sfr = np.log10(gal_sfr)

gv_mask = (gal_ssfr	> ssfr_min) & (gal_ssfr < ssfr_max)
sf_mask = gal_ssfr > ssfr_max

upper = sfms_line(sm_plot, a=1., b=ssfr_max)
lower = sfms_line(sm_plot, a=1., b=ssfr_min)

cmap = plt.get_cmap('jet_r')
new_cmap = truncate_colormap(cmap, 0.15, 0.95)

plt.plot(sm_plot, upper, ls='--', lw=1.5, c='k')
plt.plot(sm_plot, lower, ls='--', lw=1.5, c='k')
plt.axvline(10., ls='--', lw=1.5, c='k')
plt.axvline(10.5, ls='--', lw=1.5, c='k')
plt.axvline(11., ls='--', lw=1.5, c='k')
plt.scatter(gal_sm, gal_sfr, s=1, c=gal_ssfr, cmap=new_cmap)
plt.xlim(9., 11.5)
plt.clim(-12, -9)
plt.colorbar(label='log sSFR')
plt.legend()
plt.xlabel('log M*')
plt.ylabel('log SFR')
plt.savefig(plots_dir+'b18_sample_colormap.png')
plt.clf()

plt.plot(sm_plot, upper, ls='--', lw=1.5, c='k')
plt.plot(sm_plot, lower, ls='--', lw=1.5, c='k')
plt.axvline(10., ls='--', lw=1.5, c='k')
plt.axvline(10.5, ls='--', lw=1.5, c='k')
plt.axvline(11., ls='--', lw=1.5, c='k')
plt.scatter(gal_sm[gal_cent], gal_sfr[gal_cent], s=1, c=gal_ssfr[gal_cent], cmap=new_cmap)
plt.xlim(9., 11.5)
plt.clim(-12, -9)
plt.colorbar(label='log sSFR')
plt.legend()
plt.xlabel('log M*')
plt.ylabel('log SFR')
plt.savefig(plots_dir+'b18_sample_centrals_colormap.png')
plt.clf()

	

for i in range(3):
	if i != 2:
		sm_mask = (gal_sm > masses[i]) & (gal_sm < masses[i+1])
	else:
		sm_mask = gal_sm > masses[i]

	plt.plot(sm_plot, upper, ls='--', lw=1.5, c='k')
	plt.plot(sm_plot, lower, ls='--', lw=1.5, c='k')
	plt.axvline(10., ls='--', lw=1.5, c='k')
	plt.axvline(10.5, ls='--', lw=1.5, c='k')
	plt.axvline(11., ls='--', lw=1.5, c='k')
	plt.scatter(gal_sm[sm_mask*gal_cent*gv_mask], gal_sfr[sm_mask*gal_cent*gv_mask], s=1, c='g', label='Green Valley')
	plt.scatter(gal_sm[sm_mask*gal_cent*sf_mask], gal_sfr[sm_mask*gal_cent*sf_mask], s=1, c='b', label='SFMS')
	plt.xlim(9., 11.5)
	plt.legend()
	plt.xlabel('log M*')
	plt.ylabel('log SFR')
	plt.savefig(plots_dir+'b18_sample_centrals_mask_'+str(i)+'.png')
	plt.clf()

with h5py.File(h5_dir+'gv_samples.h5', 'a') as f:
	f.create_dataset(model+'_'+snap, data=np.array(gv_mask))
with h5py.File(h5_dir+'sf_samples.h5', 'a') as f:
	f.create_dataset(model+'_'+snap, data=np.array(sf_mask))
