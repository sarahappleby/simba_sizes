import matplotlib.pyplot as plt
import pygad as pg
import numpy as np
import h5py
import caesar

def magab_to_l(magab):
        """
        Find the luminosity from AB magnitude 
        """
        mab0 = 3631e-26 # conversion from janksys to W/m**2/Hz
        pc = 3.085e16 # 1 parsec in metres
        const = 4*np.pi*((pc*10)**2)*mab0
        tens = np.ones_like(magab)*10
        lum = np.power(tens, -0.4*np.array(magab))
        mask = lum != 1
        return pg.UnitArr(const*mask*lum, 'W')

model = 'm50n512'
snap = '151'
wind = 's50j7k'
plot_dir = './'
factor = 5.

gals = [45, 94, 195, 201, 259, 364]

data_dir = '/home/rad/data/'+model+'/'+wind+'/'
sim =  caesar.load(data_dir+'Groups/'+model+'_'+snap+'.hdf5')
gal_rad = np.array([i.radii['stellar_half_mass'].in_units('kpc') for i in sim.central_galaxies])

s = pg.Snap(data_dir+'snap_'+model+'_'+snap+'.hdf5')
h = s.cosmology.h()
z = s.redshift

abundance_h1 = s.gas['NeutralHydrogenAbundance']
abundance_h2 = s.gas['fh2']
mass = s.gas['mass'].in_units_of('Msol')
s.gas['mass_h1'] = abundance_h1*mass
s.gas['mass_h2'] = abundance_h2*mass

bh_pos = s.bh['pos'] / (h*(1.+z)) # in kpc

s.stars['lum_r'] = magab_to_l(np.array(s.stars['mag_r']))

for i in gals:
    pos = np.array(bh_pos[sim.central_galaxies[i].bhlist[0]])
    
    radius = str(round(factor*gal_rad[i], 2)*2) + ' kpc'
    ball_center = pg.UnitArr(pos, 'kpc')
    ball = s[pg.BallMask(radius, center=ball_center)]
    
    args = dict(cmap='magma', fontsize=8, Npx=256)
    fig, axes = plt.subplots(1,2, figsize=(9,9))
    axes = axes.flatten()
    
    pg.plotting.image(ball.gas, qty='mass_h1', xaxis=0, yaxis=1, ax=axes[0], **args)
    pg.plotting.image(ball.stars, qty='lum_r', xaxis=0, yaxis=1, ax=axes[1], **args)
    
    plt.savefig(plot_dir+'gal_'+str(i)+'_pygad.png')
    plt.clf()

