import caesar
import yt
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

def sfms_line(x0, a=0.73, b=-7.7):
	return x0*a + b

def sf_check(x0, y0, b1):
	y = sfms_line(x0, a=0.73, b=b1)
	if y < y0:
		return True
	else:
		return False

def gv_check(x0, y0, b1, b2):
	y1 = sfms_line(x0, a=0.73, b=b1)
	y2 = sfms_line(x0, a=0.73, b=b2)
	if (y1 > y0) & (y2 < y0):
		return True
	else:
		return False  

gv_ssfr_min = -13.

masses = [9., 9.5, 10., 10.5, 11.]

model = sys.argv[1]
wind = sys.argv[2]
snap = sys.argv[3]
selection = sys.argv[4]

if selection == '1':
	b1 = -7.7; b2 = -8.2
elif selection == '2':
	b1 = -7.7; b2 = -8.7
elif selection == '3':
	b1 = -7.45; b2 = -8.45

h5_dir = '/home/sapple/simba_sizes/profiles/gv_sample/belfiore/'+wind+'/selection_'+selection +'/'
if not os.path.exists(h5_dir):
	os.makedirs(h5_dir)

plots_dir = h5_dir +model+'_'+snap+'/'
if not os.path.exists(plots_dir):
	os.makedirs(plots_dir)

data_dir = '/home/rad/data/'+model+'/'+wind+'/'

sim =  caesar.load(data_dir+'Groups/'+model+'_'+snap+'.hdf5', LoadHalo=False)
c = yt.units.c.in_units('m/s')


gal_cent = np.array([i.central for i in sim.galaxies])
gal_sm = np.array([i.masses['stellar'].in_units('Msun') for i in sim.galaxies])
gal_sfr = np.array([i.sfr.in_units('Msun/yr') for i in sim.galaxies])
gal_h1_frac = np.array([i.masses['HI'].in_units('Msun') /i.masses['stellar'].in_units('Msun') for i in sim.galaxies ])
gal_h2_frac = np.array([i.masses['H2'].in_units('Msun') /i.masses['stellar'].in_units('Msun') for i in sim.galaxies ])
gal_gas_frac = np.array([i.masses['gas'].in_units('Msun') /i.masses['stellar'].in_units('Msun') for i in sim.galaxies ])
gal_mbhdot_c2 = np.array([(i.bhmdot.in_units('kg/s')*(c**2)).in_units('erg/s') for i in sim.galaxies]) * 0.1


gal_ssfr = np.log10(gal_sfr / gal_sm)
gal_sm = np.log10(gal_sm)
gal_sfr = np.log10(gal_sfr)

gal_h1_frac[np.where(gal_h1_frac == 0.)[0]] = 1.e-6
gal_h2_frac[np.where(gal_h2_frac == 0.)[0]] = 1.e-6
gal_gas_frac[np.where(gal_gas_frac == 0.)[0]] = 1.e-6

"""
Get sample of the green valley from Belfiore 18
"""
ssfr_mask = (gal_ssfr > gv_ssfr_min) 

sf_mask = [sf_check(gal_sm[i], gal_sfr[i], b1) for i in range(len(sim.galaxies))]
gv_mask = [gv_check(gal_sm[i], gal_sfr[i], b1, b2) for i in range(len(sim.galaxies))]

sm_plot = np.arange(9.0, 13., 0.5)
belfiore_main = sfms_line(sm_plot, a=0.73, b=-7.7)
if selection == '1' or selection == '2':
	belfiore_lower = sfms_line(sm_plot, a=0.73, b=b2)
elif selection == '3':
	new_upper = sfms_line(sm_plot, a=0.73, b=b1)
	new_lower = sfms_line(sm_plot, a=0.73, b=b2)

cmap = plt.get_cmap('jet_r')
new_cmap = truncate_colormap(cmap, 0.15, 0.95)

plt.plot(sm_plot, belfiore_main, ls='--', lw=1.5, c='m', label='B18')
if selection == '1':
        plt.plot(sm_plot, belfiore_lower, ls='-.', lw=1.5, c='m', label='B18 - 0.5 dex')
elif selection == '2':
        plt.plot(sm_plot, belfiore_lower, ls='--', lw=1.5, c='m')
else:
        plt.plot(sm_plot, new_upper, ls='--', lw=1.5, c='b', label='B18 + 0.25 dex')
        plt.plot(sm_plot, new_lower, ls='-.', lw=1.5, c='b', label='B18 - 0.75 dex')
plt.axvline(10., ls=':', lw=1.8, c='k')
plt.axvline(10.5, ls=':', lw=1.8, c='k')
plt.axvline(11., ls=':', lw=1.8, c='k')
plt.annotate('Low', xy=(0.3, 0.9), xycoords='axes fraction',size=16)
plt.annotate('Int', xy=(0.46, 0.9), xycoords='axes fraction',size=16)
plt.annotate('High', xy=(0.7, 0.9), xycoords='axes fraction',size=16)
plt.scatter(gal_sm, gal_sfr, s=0.3, c=np.log10(gal_mbhdot_c2), cmap=new_cmap)
plt.xlim(9.0,12.5)
plt.ylim(-3.5, )
plt.clim(40, 44)
plt.colorbar(label=r'$ \textrm{log}\ (0.1 \dot{M}_{bh}c^2 / \textrm{erg s}^{-1})$')
plt.legend(fontsize=12, loc=4)
plt.xlabel(r'$\log\ (M_{*} / M_{\odot})$')
plt.ylabel(r'$\textrm{log} (\textrm{SFR} / M_{\odot}\textrm{yr}^{-1})$')
plt.savefig(plots_dir+'b18_sample_all_colormap_mbhdot.png')
plt.clf()



plt.plot(sm_plot, belfiore_main, ls='--', lw=1.5, c='m', label='B18')
if selection == '1':
        plt.plot(sm_plot, belfiore_lower, ls='-.', lw=1.5, c='m', label='B18 - 0.5 dex')
elif selection == '2':
        plt.plot(sm_plot, belfiore_lower, ls='--', lw=1.5, c='m')
else:
        plt.plot(sm_plot, new_upper, ls='--', lw=1.5, c='b', label='B18 + 0.25 dex')
        plt.plot(sm_plot, new_lower, ls='-.', lw=1.5, c='b', label='B18 - 0.75 dex')
plt.axvline(10., ls=':', lw=1.8, c='k')
plt.axvline(10.5, ls=':', lw=1.8, c='k')
plt.axvline(11., ls=':', lw=1.8, c='k')
plt.annotate('Low', xy=(0.3, 0.9), xycoords='axes fraction',size=16)
plt.annotate('Int', xy=(0.46, 0.9), xycoords='axes fraction',size=16)
plt.annotate('High', xy=(0.7, 0.9), xycoords='axes fraction',size=16)
plt.scatter(gal_sm, gal_sfr, s=0.3, c=np.log10(gal_h1_frac), cmap=new_cmap)
plt.xlim(9.0,12.5)
plt.ylim(-3.5, )
#plt.clim(-12, -9)
#plt.colorbar(label=r'$\textrm{log} (\textrm{sSFR} / \textrm{yr}^{-1})$')
plt.clim(-2, 0.5)
plt.colorbar(label=r'$ \textrm{log}\ f_{HI}$')
plt.legend(fontsize=12, loc=4)
plt.xlabel(r'$\log\ (M_{*} / M_{\odot})$')
plt.ylabel(r'$\textrm{log} (\textrm{SFR} / M_{\odot}\textrm{yr}^{-1})$')
plt.savefig(plots_dir+'b18_sample_all_colormap_fh1.png')
plt.clf()


plt.plot(sm_plot, belfiore_main, ls='--', lw=1.5, c='m', label='B18')
if selection == '1':
        plt.plot(sm_plot, belfiore_lower, ls='-.', lw=1.5, c='m', label='B18 - 0.5 dex')
elif selection == '2':
	plt.plot(sm_plot, belfiore_lower, ls='--', lw=1.5, c='m')
else:
	plt.plot(sm_plot, new_upper, ls='--', lw=1.5, c='b', label='B18 + 0.25 dex')
	plt.plot(sm_plot, new_lower, ls='-.', lw=1.5, c='b', label='B18 - 0.75 dex')
plt.axvline(10., ls=':', lw=1.8, c='k')
plt.axvline(10.5, ls=':', lw=1.8, c='k')
plt.axvline(11., ls=':', lw=1.8, c='k')
plt.annotate('Low', xy=(0.3, 0.9), xycoords='axes fraction',size=16)
plt.annotate('Int', xy=(0.46, 0.9), xycoords='axes fraction',size=16)
plt.annotate('High', xy=(0.7, 0.9), xycoords='axes fraction',size=16)
plt.scatter(gal_sm, gal_sfr, s=0.3, c=np.log10(gal_h2_frac), cmap=new_cmap)
plt.xlim(9.0,12.5)
plt.ylim(-3.5, )
#plt.clim(-12, -9)
#plt.colorbar(label=r'$\textrm{log} (\textrm{sSFR} / \textrm{yr}^{-1})$')
plt.clim(-2, 0.5)
plt.colorbar(label=r'$ \textrm{log}\ f_{H_2}$')
plt.legend(fontsize=12, loc=4)
plt.xlabel(r'$\log\ (M_{*} / M_{\odot})$')
plt.ylabel(r'$\textrm{log} (\textrm{SFR} / M_{\odot}\textrm{yr}^{-1})$')
plt.savefig(plots_dir+'b18_sample_all_colormap_fh2.png')
plt.clf()

plt.plot(sm_plot, belfiore_main, ls='--', lw=1.5, c='m', label='B18')
if selection == '1':
        plt.plot(sm_plot, belfiore_lower, ls='-.', lw=1.5, c='m', label='B18 - 0.5 dex')
elif selection == '2':
        plt.plot(sm_plot, belfiore_lower, ls='--', lw=1.5, c='m')
else:
        plt.plot(sm_plot, new_upper, ls='--', lw=1.5, c='b', label='B18 + 0.25 dex')
        plt.plot(sm_plot, new_lower, ls='-.', lw=1.5, c='b', label='B18 - 0.75 dex')
plt.axvline(10., ls=':', lw=1.8, c='k')
plt.axvline(10.5, ls=':', lw=1.8, c='k')
plt.axvline(11., ls=':', lw=1.8, c='k')
plt.annotate('Low', xy=(0.3, 0.9), xycoords='axes fraction',size=16)
plt.annotate('Int', xy=(0.46, 0.9), xycoords='axes fraction',size=16)
plt.annotate('High', xy=(0.7, 0.9), xycoords='axes fraction',size=16)
plt.scatter(gal_sm, gal_sfr, s=0.3, c=gal_ssfr, cmap=new_cmap)
plt.xlim(9.0,12.5)
plt.ylim(-3.5, )
plt.clim(-12, -9)
plt.colorbar(label=r'$\textrm{log} (\textrm{sSFR} / \textrm{yr}^{-1})$')
plt.legend(fontsize=12, loc=4)
plt.xlabel(r'$\log\ (M_{*} / M_{\odot})$')
plt.ylabel(r'$\textrm{log} (\textrm{SFR} / M_{\odot}\textrm{yr}^{-1})$')
plt.savefig(plots_dir+'b18_sample_all_colormap_ssfr.png')
plt.clf()


with h5py.File(h5_dir+'gv_samples.h5', 'a') as f:
	f.create_dataset(model+'_'+snap, data=np.array(gv_mask))
with h5py.File(h5_dir+'sf_samples.h5', 'a') as f:
	f.create_dataset(model+'_'+snap, data=np.array(sf_mask))
