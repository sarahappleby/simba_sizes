import h5py
import numpy as np
import matplotlib.pyplot as plt
from plotting_methods import plot_all
import sys

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 14})

model = sys.argv[1]
wind = sys.argv[2]
snap = sys.argv[3]

basic_dir = '/home/sapple/simba_sizes/profiles/ssfr/extended_profiles/'
gv_centrals_dir = basic_dir + 'centrals/'+model+'_'+snap+'/'+wind+'/green_valley/random_orientation/'
gv_sats_dir = basic_dir + 'satellites/'+model+'_'+snap+'/'+wind+'/green_valley/random_orientation/'
sf_centrals_dir = basic_dir + 'centrals/'+model+'_'+snap+'/'+wind+'/star_forming/random_orientation/'
sf_sats_dir = basic_dir + 'satellites/'+model+'_'+snap+'/'+wind+'/star_forming/random_orientation/'
data_dirs = [sf_centrals_dir, gv_centrals_dir, sf_sats_dir, gv_sats_dir]

results_dir = '/home/sapple/simba_sizes/profiles/plotting/'


xlim = [0., 5.]
ylim = None
xlabel = r'$R_{half}$'

# h1:

filename = 'h1_frac_data.h5'
ylabel = r'$\textrm{log} (\Sigma_{f\ HI} / \textrm{kpc}^{-2})$'
savefile = results_dir+'h1_frac_'+wind+'.png'
plot_all(data_dirs, filename, xlabel, ylabel, xlim, ylim, savefile)

filename = 'h1_mass_data.h5'
ylabel = r'$ \textrm{log} (\Sigma_{HI} / M_{\odot}\textrm{kpc}^{-2})$'
savefile = results_dir+'h1_mass_'+wind+'.png'
plot_all(data_dirs, filename, xlabel, ylabel, xlim, ylim, savefile)

filename = 'h1_mass_sm_data.h5'
ylabel = r'$\textrm{log} (M_{HI} / M_*)$'
savefile = results_dir+'h1_mass_sm_'+wind+'.png'
plot_all(data_dirs, filename, xlabel, ylabel, xlim, ylim, savefile)

# h2:

filename = 'h2_frac_data.h5'
ylabel = r'$\textrm{log} (\Sigma_{f\ H_2} / \textrm{kpc}^{-2})$'
savefile = results_dir+'h2_frac_'+wind+'.png'
plot_all(data_dirs, filename, xlabel, ylabel, xlim, ylim, savefile)

filename = 'h2_mass_data.h5'
ylabel = r'$\textrm{log} (\Sigma_{H_2} / M_{\odot}\textrm{kpc}^{-2})$'
savefile = results_dir+'h2_mass_'+wind+'.png'
plot_all(data_dirs, filename, xlabel, ylabel, xlim, ylim, savefile)

filename = 'h2_mass_sm_data.h5'
ylabel = r'$\textrm{log} (M_{H_2} / M_*)$'
savefile = results_dir+'h2_mass_sm_'+wind+'.png'
plot_all(data_dirs, filename, xlabel, ylabel, xlim, ylim, savefile)


