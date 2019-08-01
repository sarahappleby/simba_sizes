import matplotlib.pyplot as plt
import pygad as pg
import numpy as np
import h5py
import caesar

ssfr_lim = -1.5
n_lim = 1000

model = 'm50n512'
snap = '151'
wind = 's50j7k'
plot_dir = './all_pygad_plots/'
factor = 10.
softening = pg.UnitArr([0., 0., 0., 0., 0.25, 0.], 'ckpc h_0**-1')

xaxis = [0, 1, 2]
yaxis = [1, 2, 0]
zaxis = [2, 0, 1]

gals = [20, 138, 230]

data_dir = '/home/rad/data/'+model+'/'+wind+'/'
sim =  caesar.load(data_dir+'Groups/'+model+'_'+snap+'.hdf5')

# stellar mass maps
gal_n = np.array([len(i.glist) for i in sim.central_galaxies])
gal_rad = np.array([i.radii['stellar_half_mass'].in_units('kpc') for i in sim.central_galaxies])
gal_sm = np.array([i.masses['stellar'].in_units('Msun') for i in sim.central_galaxies])
gal_sfr = np.array([i.sfr.in_units('Msun/Gyr') for i in sim.central_galaxies])
gal_ssfr = np.log10(gal_sfr / gal_sm)

ssfr_mask = gal_ssfr > ssfr_lim
n_mask = gal_n > n_lim
gals = np.arange(len(sim.central_galaxies))[(ssfr_mask*n_mask)]

s = pg.Snap(data_dir+'snap_'+model+'_'+snap+'.hdf5')
h = s.cosmology.h()
z = s.redshift

bh_pos = s.bh['pos'].in_units_of('kpc')
abundance_h1 = s.gas['NeutralHydrogenAbundance']
abundance_h2 = s.gas['fh2']
mass = s.gas['mass'].in_units_of('Msol')
s.gas['mass_h1'] = abundance_h1*mass
s.gas['mass_h2'] = abundance_h2*mass

for i in gals:
    pos = bh_pos[sim.central_galaxies[i].bhlist[0]]

    radius = str(round((factor+2)*gal_rad[i], 2)) + ' kpc'
    ball = s[pg.BallMask(radius, center=pos)]
   
    #extent = pg.UnitArr([[pos[0] - factor*gal_rad[i], pos[0] + factor*gal_rad[i]], [pos[1] - factor*gal_rad[i], pos[1] + factor*gal_rad[i]]], 'kpc')
    #args_sm = dict(cmap='magma', fontsize=8, Npx=256, softening=softening, extent=extent, vlim=[10.**-3.5, 10.**-0.1])
    #args_h1 = dict(cmap='magma', fontsize=8, Npx=256, extent=extent)

    args_sm = dict(cmap='magma', fontsize=8, Npx=1024, softening=softening, vlim=[10.**-3.5, 10.**-0.1])
    args_h1 = dict(cmap='magma', fontsize=8, Npx=1024, vlim=[4.e3, 1.e9], cbartitle='')

    fig, axes = plt.subplots(1,3, figsize=(12,6))
    for j in range(len(xaxis)):
        pg.plotting.image(ball.gas, qty='mass_h1', xaxis=xaxis[j], yaxis=yaxis[j], ax=axes[j], **args_h1)
    plt.savefig(plot_dir+'gal_'+str(i)+'_h1_pygad.png')
    plt.clf()
   
    """
    radius = str(round((factor)*gal_rad[i], 2)) + ' kpc'
    ball = s[pg.BallMask(radius, center=pos)]
    args_h1 = dict(cmap='magma', fontsize=8, Npx=1024, vlim=[10.**3., 10.**9.])

    fig, axes = plt.subplots(1,3, figsize=(12,6))
    for j in range(len(xaxis)):
        pg.plotting.image(ball.gas, qty='mass_h2', xaxis=xaxis[j], yaxis=yaxis[j], ax=axes[j], **args_h1)
    plt.savefig(plot_dir+'gal_'+str(i)+'_h2_small_pygad.png')
    plt.clf()
    """
"""
# HI column density
gal_rad = np.array([i.radii['stellar_half_mass'].in_units('cm') for i in sim.central_galaxies])

s = pg.Snap(data_dir+'snap_'+model+'_'+snap+'.hdf5')
h = s.cosmology.h()
z = s.redshift

abundance_h1 = s.gas['NeutralHydrogenAbundance']
mass = s.gas['mass'].in_units_of('Msol')
s.gas['mass_h1'] = abundance_h1*mass
s.gas['N_HI'] = s.gas['mass_h1'].in_units_of('1e+30 kg') / pg.cosmology.m_p.in_units_of('kg')

bh_pos = s.bh['pos'].in_units_of('cm')

s['pos'].convert_to('cm')

for i in gals:
    pos = bh_pos[sim.central_galaxies[i].bhlist[0]]
    
    radius = str(round(factor*gal_rad[i], 2)) + ' cm'
    ball = s[pg.BallMask(radius, center=pos)]
   
    extent = pg.UnitQty(gal_rad[i]*factor, 'cm')
    extent = [[pg.UnitQty(pos[0] - factor*gal_rad[i], 'cm'), pg.UnitQty(pos[0] + factor*gal_rad[i], 'cm')], 
             [pg.UnitQty(pos[1] - factor*gal_rad[i], 'cm'), pg.UnitQty(pos[1] + factor*gal_rad[i], 'cm')]]

    args = dict(cmap='magma', fontsize=8, Npx=256)
    ball.gas['pos'].convert_to('cm')
    ball.gas['hsml'].convert_to('cm')
    pg.plotting.image(ball.gas, qty='N_HI', xaxis=0, yaxis=1, **args)
    plt.savefig(plot_dir+'gal_'+str(i)+'_HI_pygad.png')
    plt.clf()
"""
