import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

import numpy as np

import caesar

import os
import sys
sys.path.append('/home/sapple/tools/')
from subdivide import *

data_dir = '/home/sapple/simba_sizes/data/'
plot_dir = '/home/sapple/simba_sizes/plots/'
mufasa_caesar = '/home/rad/data/m50n512/fh_qr/Groups/'
simba_caesar = '/home/rad/data/m50n512/s50j7k/Groups/'

m_star = 5.e10
m_star_str = str(m_star)[0] + 'e' + str(int(np.log10(m_star)))

mufasa_snaps = ['070', '085', '095', '105', '125', '126', '135']
simba_snaps = ['062', '078', '090', '105', '125', '145', '151']

mufasa_z = [3.0, 2.0, 1.5, 1.0, 0.5, 0.2, 0.0]
simba_z = [3.0, 2.0, 1.5, 1.0, 0.5, 0.1, 0.0]

mufasa_sizes = np.zeros(len(mufasa_z))
mufasa_per_25 = np.zeros(len(mufasa_z))
mufasa_per_75 = np.zeros(len(mufasa_z))
simba_sizes = np.zeros(len(mufasa_z))
simba_per_25 = np.zeros(len(mufasa_z))
simba_per_75 = np.zeros(len(mufasa_z))

i = 0
for snap in mufasa_snaps:
	if not os.path.isfile(data_dir+'s_mass_m50n512_fh_qr_'+snap+'.npy'):
		sim = caesar.load(mufasa_caesar+'m50n512_'+snap+'.hdf5')
		s_mass = np.array([j.masses['stellar'] for j in sim.galaxies])
		np.save(data_dir+'s_mass_m50n512_fh_qr_'+snap+'.npy', s_mass)
	else:
		s_mass = np.load(data_dir+'s_mass_m50n512_fh_qr_'+snap+'.npy')
	mask = (s_mass > 0.95*m_star) & (s_mass < 1.05*m_star)
	r = np.load(data_dir+'halflight_m50n512_fh_qr_'+snap+'.npy')
	mufasa_sizes[i] = np.median(r); mufasa_per_25[i] = np.percentile(r, 25); mufasa_per_75[i] = np.percentile(r, 75)
	i += 1

i = 0
for snap in simba_snaps:
        if not os.path.isfile(data_dir+'s_mass_m50n512_s50j7k_'+snap+'.npy'):
                sim = caesar.load(simba_caesar+'m50n512_'+snap+'.hdf5')
                s_mass = np.array([j.masses['stellar'] for j in sim.galaxies])
                np.save(data_dir+'s_mass_m50n512_s50j7k_'+snap+'.npy', s_mass)
        else:
                s_mass = np.load(data_dir+'s_mass_m50n512_s50j7k_'+snap+'.npy')
        mask = (s_mass > 0.95*m_star) & (s_mass < 1.05*m_star)
        r = np.load(data_dir+'halflight_m50n512_s50j7k_'+snap+'.npy')
        simba_sizes[i] = np.median(r); simba_per_25[i] = np.percentile(r, 25); simba_per_75[i] = np.percentile(r, 75)
        i += 1	

plt.plot(mufasa_z, mufasa_sizes, linestyle='-', c='g', label='Mufasa')
plt.fill_between(mufasa_z, mufasa_per_25, mufasa_per_75, facecolor='g', alpha=0.15, linewidth=1)

plt.plot(simba_z, simba_sizes, linestyle='-', c='b', label='Simba')
plt.fill_between(simba_z, simba_per_25, simba_per_75, facecolor='b', alpha=0.15, linewidth=1)

plt.xlabel(r'z', fontsize=16)
plt.xlim(3.1, -0.1)
plt.ylabel(r'$\log\ R_{half,*}$' ,fontsize=16)
plt.legend(loc=2)

plt.savefig(plot_dir+'redshift_size_'+m_star_str+'.png')
plt.clf()
