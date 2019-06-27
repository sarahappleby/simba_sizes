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

xlim = [0., 3.]
ylim = None
xlabel = r'$R_{half}$'

filename = 'ages_data.h5'
ylabel = r'$ \textrm{log} (\tau_{*} / \textrm{Myr})$'
savefile = results_dir+'ages_'+wind+'.png'
plot_all(data_dirs, filename, xlabel, ylabel, xlim, ylim, savefile)
