import caesar
from readgadget import *
import sys
import numpy as np
import h5py
# ---------------------------
def get_r(mass,pos,r):
    # positions need to be in kpc
    mtot = np.sum(mass)
    pos_new = np.zeros_like(pos)

    x0 = np.sum(pos[:,0]*mass)/np.sum(mass)
    y0 = np.sum(pos[:,1]*mass)/np.sum(mass)
    z0 = np.sum(pos[:,2]*mass)/np.sum(mass)
    pos_new[:,0] = pos[:,0] - x0
    pos_new[:,1] = pos[:,1] - y0
    pos_new[:,2] = pos[:,2] - z0

    mr50 = 0.
    dr = 0.1              # in kpc
    i = 0
    while (mr50 < mtot*r/100.):
         i += 1
         mask = ((pos_new[:,0]**2 + pos_new[:,1]**2 + pos_new[:,2]**2) <= (dr*i)**2)
         mr50 = np.sum(mass[mask])
    return dr*i,x0,y0,z0

def get_cm(mass,pos):
    # positions need to be in kpc
    mtot = np.sum(mass)
    pos_new = np.zeros_like(pos)

    x0 = np.sum(pos[:,0]*mass)/np.sum(mass)
    y0 = np.sum(pos[:,1]*mass)/np.sum(mass)
    z0 = np.sum(pos[:,2]*mass)/np.sum(mass)

    return x0,y0,z0


def get_cm_radius(pos,mass,dx):
    #d_one = 30.  # in kpc
    one_galaxy = 0
    dr = 5.        # step in radius in kpc
    dens_threshold = 5   # if the first max is dens_threshold-time denser than the second one, there is only one galaxy 

    x_min = np.min(pos[:,0])
    x_max = np.max(pos[:,0])
    y_min = np.min(pos[:,1])
    y_max = np.max(pos[:,1])
    z_min = np.min(pos[:,2])
    z_max = np.max(pos[:,2])

    delta_x = x_max-x_min
    delta_y = y_max-y_min
    delta_z = z_max-z_min
    d_min = 0.5*np.min([delta_x,delta_y,delta_z])
    d_min2 = 0.5*d_min
    print 'd_min,d_min2:',d_min,d_min2

    binx = np.int((x_max - x_min)/dx)
    biny = np.int((y_max - y_min)/dx)
    binz = np.int((z_max - z_min)/dx)

    #binx = 50
    #biny = binx
    #binz = binx
    h, [x,y,z] = np.histogramdd((pos[:,0],pos[:,1],pos[:,2]),bins=(binx,biny,binz),weights=mass)

    i,j,k = np.unravel_index(h.argmax(), h.shape)
    h_max1 = h[i,j,k]
    x_max1 = x[i]
    y_max1 = y[j]
    z_max1 = z[k]

    # cm1 corrected
    d_min_i = 5  # initial guess for radius in kpc
    stop = 0
    kk = 0
    cm_x_i = x_max1
    cm_y_i = y_max1
    cm_z_i = z_max1

    while (stop == 0):
         mask1 = (np.sqrt((posstar[:,0] - cm_x_i)**2+(posstar[:,1] - cm_y_i)**2+(posstar[:,2] - cm_z_i)**2) <= d_min_i*(kk + 1))
         mass_temp1 = np.sum(mstar[mask1])
         cm_x_olds1 = np.sum(posstar[:,0][mask1]*mstar[mask1])/mass_temp1
         cm_y_olds1 = np.sum(posstar[:,1][mask1]*mstar[mask1])/mass_temp1
         cm_z_olds1 = np.sum(posstar[:,2][mask1]*mstar[mask1])/mass_temp1
        
         if ( ((np.abs(cm_x_olds1 - cm_x_i)) < 0.01*cm_x_olds1) and ((np.abs(cm_y_olds1 - cm_y_i)) < 0.01*cm_y_olds1) and ((np.abs(cm_z_olds1 - cm_z_i)) < 0.01*cm_z_olds1) ):
           stop = 1
         else:
           cm_x_i = cm_x_olds1
           cm_y_i = cm_y_olds1
           cm_z_i = cm_z_olds1
           kk += 1
    print 'cm:',cm_x_olds1,cm_y_olds1,cm_z_olds1
    ######################################################
    # mask out the galaxy 1
    ii_x = int(d_min*(1 + kk)/((x_max-x_min)/binx))
    ii_y = int(d_min*(1 + kk)/((y_max-y_min)/biny))
    ii_z = int(d_min*(1 + kk)/((z_max-z_min)/binz))

    ii_max = np.max([ii_x,ii_y,ii_z])
    ii_one_x = int(d_min2/((x_max-x_min)/binx))
    ii_one_y = int(d_min2/((y_max-y_min)/biny))
    ii_one_z = int(d_min2/((z_max-z_min)/binz))

    ii_one_min = np.min([ii_one_x,ii_one_y,ii_one_z])
    ii = np.min([ii_max,ii_one_min])
    h[i-ii:i+ii,j-ii:j+ii,k-ii:k+ii] = 0.0

    i,j,k = np.unravel_index(h.argmax(), h.shape)
    h_max2 = h[i,j,k]
    x_max2 = x[i]
    y_max2 = y[j]
    z_max2 = z[k]

    print 'h_max1,h_max2:',h_max1,h_max2
    print 'dist:',np.sqrt((x_max2-x_max1)**2+(y_max2-y_max1)**2+(z_max2-z_max1)**2)

    if (h_max1 > dens_threshold*h_max2 or np.sqrt((x_max2-x_max1)**2+(y_max2-y_max1)**2+(z_max2-z_max1)**2) <= d_min2):
      one_galaxy=1       # there is only 1 galaxy
      print 'one galaxy'
      """mask = (np.sqrt((pos[:,0]-x_max1)**2+(pos[:,1]-y_max1)**2+(pos[:,2]-z_max1)**2)<=(d_min))
      mass_temp = np.sum(mstar[mask])
      cm_x_olds1 = np.sum(pos[:,0][mask]*mstar[mask])/mass_temp
      cm_y_olds1 = np.sum(pos[:,1][mask]*mstar[mask])/mass_temp
      cm_z_olds1 = np.sum(pos[:,2][mask]*mstar[mask])/mass_temp
      """
    else:
      #########################
      # cm2 corrected
      #########################
      d_min2_i = d_min/2.  # initial guess for radius in kpc
      stop = 0
      kk2 = 0
      cm_x_i = x_max2
      cm_y_i = y_max2
      cm_z_i = z_max2
      print 'cm2:',cm_x_i,cm_y_i,cm_z_i

      while (stop == 0):
         mask1 = (np.sqrt((pos[:,0] - cm_x_i)**2+(pos[:,1] - cm_y_i)**2+(pos[:,2] - cm_z_i)**2) <= d_min2_i*(kk2 + 1))
         if (len(mstar[mask1]) > 0):
            mass_temp1 = np.sum(mstar[mask1])
            cm_x_olds2 = np.sum(pos[:,0][mask1]*mstar[mask1])/mass_temp1
            cm_y_olds2 = np.sum(pos[:,1][mask1]*mstar[mask1])/mass_temp1
            cm_z_olds2 = np.sum(pos[:,2][mask1]*mstar[mask1])/mass_temp1
        
            if ( ((np.abs(cm_x_olds2 - cm_x_i)) < 0.01*cm_x_olds2) and ((np.abs(cm_y_olds2 - cm_y_i)) < 0.01*cm_y_olds2) and ((np.abs(cm_z_olds2 - cm_z_i)) < 0.01*cm_z_olds2) ):
              stop = 1
            else:
              cm_x_i = cm_x_olds2
              cm_y_i = cm_y_olds2
              cm_z_i = cm_z_olds2
            kk2 += 1
         else:
           cm_x_olds2 = cm_x_i
           cm_y_olds2 = cm_y_i
           cm_z_olds2 = cm_z_i
           stop = 1

      #########################
      d_maxs = np.sqrt((cm_x_olds2-cm_x_olds1)**2+(cm_y_olds2-cm_y_olds1)**2+(cm_z_olds2-cm_z_olds1)**2)
      d_maxs_2D = np.sqrt((cm_x_olds2-cm_x_olds1)**2+(cm_y_olds2-cm_y_olds1)**2)
      print 'd_maxs:',d_maxs

      ##############################
      # half-mass radius2 (olds) 3D
      ##############################
      if (d_maxs <= d_min2 or dr >= d_maxs):
        one_galaxy = 1
      else:
        dd = 0
        dens_vol = []
     
        while ((dr*(1 + dd)) < d_maxs): 
          mask2 = (np.sqrt((pos[:,0]-cm_x_olds2)**2+(pos[:,1]-cm_y_olds2)**2+(pos[:,2]-cm_z_olds2)**2) <= dr*(1 + dd))
          mass_temp = np.sum(mstar[mask2] )
          if (dd == 0):
            volume_i = 4/3.*np.pi*dr**3
          else:
            volume_i = 4/3.*np.pi*dr**3*((1+dd)**3 - dd**3 )

          dens_i = mass_temp #/volume_i  # don't use volume i.e. 3D density (as it decreases rapidly because of oncreasing volume), mass seems better
          dens_vol.append(dens_i)
          dd += 1

        dens_vol = np.asarray(dens_vol)
        dens_vol /= 10**9

        diff = np.zeros(len(dens_vol) - 1)
        arr1 = dens_vol[1:]
        arr2 = dens_vol[:len(dens_vol)-1]
        diff = (arr1 - arr2)/arr1*100.
        frac_diff = 1.
        jj = 0
        idx = []
      
        if (len(arr1) > 0):
          while (len(idx) == 0):
             #print 'jj:',jj,frac_diff*(1+jj)
             idx = np.where(diff < frac_diff*(1+jj))[0]
             jj += 1

          radius2 = dr*(idx[0])
        else:
          one_galaxy = 1
          radius2 = -99

        print 'radius2:',radius2
    ##############################
    # half-mass radius1 (olds) 3D
    ##############################
    # only if there are 2 galaxies
    if (one_galaxy == 0):
        dd = 0
        dens_vol = []

        while ((dr*(1 + dd)) < d_maxs):
            mask2 = (np.sqrt((pos[:,0]-cm_x_olds1)**2+(pos[:,1]-cm_y_olds1)**2+(pos[:,2]-cm_z_olds1)**2) <= dr*(1 + dd))
            mass_temp = np.sum(mstar[mask2] )
            if (dd == 0):
              volume_i = 4/3.*np.pi*dr**3
            else:
              volume_i = 4/3.*np.pi*dr**3*((1+dd)**3 - dd**3 )

            dens_i = mass_temp #/volume_i
            dens_vol.append(dens_i)
            dd += 1

        dens_vol = np.asarray(dens_vol)
        dens_vol /= 10**9

        diff = np.zeros(len(dens_vol) - 1)
        arr1 = dens_vol[1:]
        arr2 = dens_vol[:len(dens_vol)-1]
        diff = (arr1 - arr2)/arr1*100.

        frac_diff = 1.
        jj = 0
        idx = []
        while (len(idx) == 0):
           idx = np.where(diff < frac_diff*(1+jj))[0]
           jj += 1
        radius = dr*(idx[0])

    if (one_galaxy == 1):
      radius = -99

    print 'radius:',radius

    if (one_galaxy == 1):
      return cm_x_olds1,cm_y_olds1,cm_z_olds1,radius
    else:
      return cm_x_olds1,cm_y_olds1,cm_z_olds1,radius, cm_x_olds2,cm_y_olds2,cm_z_olds2,radius2


def get_ang_mtm(mask, pos, vel, mass):

    galPosx = np.sum(mstar[mask]*posstar[:,0][mask])/np.sum(mstar[mask])
    galPosy = np.sum(mstar[mask]*posstar[:,1][mask])/np.sum(mstar[mask])
    galPosz = np.sum(mstar[mask]*posstar[:,2][mask])/np.sum(mstar[mask])

    print 'galPos:',galPosx,galPosy,galPosz

    galVelx = np.sum(mstar[mask]*velstar[:,0][mask])/np.sum(mstar[mask])
    galVely = np.sum(mstar[mask]*velstar[:,1][mask])/np.sum(mstar[mask])
    galVelz = np.sum(mstar[mask]*velstar[:,2][mask])/np.sum(mstar[mask])

    vx = velstar[:,0][mask] - galVelx
    vy = velstar[:,1][mask] - galVely
    vz = velstar[:,2][mask] - galVelz

    x = posstar[:,0][mask] - galPosx #* kpc2km    # in km; not really needed <=> normalised afterwards
    y = posstar[:,1][mask] - galPosy #* kpc2km
    z = posstar[:,2][mask] - galPosz #* kpc2km

    Lx = np.sum(mstar[mask]*(y*vz - z*vy))   # in Msun * km^{2} * s^{-1}
    Ly = np.sum(mstar[mask]*(z*vx - x*vz))
    Lz = np.sum(mstar[mask]*(x*vy - y*vx))

    Lnorm = np.sqrt(Lx**2 + Ly**2 + Lz**2)

    if (Lnorm > 0):
        Lx /= Lnorm
        Ly /= Lnorm
        Lz /= Lnorm

        # cylindrical coordinates:
        l = x*Lx + y*Ly + z*Lz    #x,y,z not posstar 

        L0_x = galPosx + l*Lx
        L0_y = galPosy + l*Ly
        L0_z = galPosz + l*Lz

        # ez = (Lxgal,Lygal,Lzgal)
        # er = (rx,ry,rz)
        # etheta = (thetax,thetay,thetaz) = ez x er

        # unit radial vector to particles
        rx = posstar[:,0][mask] - L0_x
        ry = posstar[:,1][mask] - L0_y
        rz = posstar[:,2][mask] - L0_z
        rnorm = np.sqrt(rx**2 + ry**2 +rz**2)

        rx[rnorm > 0] /= rnorm[rnorm > 0]
        ry[rnorm > 0] /= rnorm[rnorm > 0]
        rz[rnorm > 0] /= rnorm[rnorm > 0]

        thetax = Ly*rz - Lz*ry
        thetay = Lz*rx - Lx*rz 
        thetaz = Lx*ry - Ly*rx

        # velocity components in the cylindrical coords
        vcz = vx*Lx + vy*Ly + vz*Lz
        vcr = vx*rx + vy*ry + vz*rz
        vct = vx*thetax + vy*thetay + vz*thetaz

        vrmean0 = np.sum(vcr*mstar[mask])/np.sum(mstar[mask])
        vzmean0 = np.sum(vcz*mstar[mask])/np.sum(mstar[mask])
        vtmean0 = np.sum(vct*mstar[mask])/np.sum(mstar[mask])

        sigrmean0 = np.sqrt(np.sum(( vcr - vrmean0)**2*mstar[mask])/np.sum(mstar[mask]))
        sigzmean0 = np.sqrt(np.sum(( vcz - vzmean0)**2*mstar[mask])/np.sum(mstar[mask]))
        sigtmean0 = np.sqrt(np.sum(( vct - vtmean0)**2*mstar[mask])/np.sum(mstar[mask]))

    else:
        Lx, Ly, Lz, vrmean0, vzmean0, vtmean0, sigrmean0, sigzmean0, sigtmean0 = [-99]*9

    return Lx, Ly, Lz, vrmean0, vzmean0, vtmean0, sigrmean0, sigzmean0, sigtmean0
        

# ---------------------------
snap = '151'   # z=0: 151; z=0.5: 125; z=1: 105; z=2: 078
model = 'm100n1024'
output = './Gals_'+model+'_'+str(snap)+'.asc'

# simba
snapfile = '/home/rad/data/'+model+'/s50j7k/snap_'+model+'_'+str(snap)+'.hdf5'
caesarfile = '/home/rad/data/'+model+'/s50j7k/Groups/'+model+'_'+str(snap)+'.hdf5'

# load in input file
sim = caesar.load(caesarfile,LoadHalo=False)  # load without halos
ngal = sim.ngalaxies
print 'ngal:',ngal

#boxsize
bsize  = sim.simulation.boxsize.to('kpccm').d
print 'boxsize',bsize

redshift = np.round(sim.simulation.redshift,decimals=2)
print 'redshift:',redshift
h = sim.simulation.hubble_constant
print 'h:',h
kpc2km = 3.086*10**(16.)     # conversion from kpc to km


# read in the particle information we need, mass, positions and velocity
"""pmass = readsnap(snapfile,'mass','gas',units=1,suppress=1)/h  # note the h^-1 from Gadget units
ppos = readsnap(snapfile,'pos','gas',units=0,suppress=1)/h      #/h/1000.-> in Mpc; otherwise in kpc comoving: max 50Mpc/h
pvel = readsnap(snapfile,'vel','gas',units=0,suppress=1)  # in km/s
h1 = readsnap(snapfile,'NeutralHydrogenAbundance','gas',units=0,suppress=1)  # in km/s
"""
gaspos = readsnap(snapfile,'pos','gas',units=0,suppress=1)/h
pmass = readsnap(snapfile,'mass','star',units=1,suppress=1)/h  # note the h^-1 from Gadget units
ppos = readsnap(snapfile,'pos','star',units=0,suppress=1)/h      #/h/1000.-> in Mpc; otherwise in kpc comoving: max 50Mpc/h
pvel = readsnap(snapfile,'vel','star',units=0,suppress=1)  # in km/s

print 'min-max ppos:',np.min(ppos),np.max(ppos)
print 'min-max pmass:',np.min(pmass),np.max(pmass)

# compile the galaxy data from the caesar file
gals = np.asarray([i for i in sim.galaxies])
gmass = np.asarray([i.masses['stellar'] for i in sim.galaxies])
gpos = np.asarray([i.pos for i in sim.galaxies])    # kpc comoving
gvel = np.asarray([i.vel for i in sim.galaxies])

galID = np.asarray([i.GroupID for i in sim.galaxies]).astype(np.int)

print 'min-max gpos:',np.min(gpos),np.max(gpos)

# compute the angular momentum from the particle lists
Lxgal = []
Lygal = []
Lzgal = []

vrmean = []
vzmean = []
vtmean = []
mtot = []
sigrmean = []
sigzmean = []
sigtmean = []

r50_tot = []
nb_part = []

slists = []
glists = []
new_pos = []
new_rad = []
gal_ids = []

t = 0
dx = 1.  # bin size in [kpc] to get the max density

"""gmass = np.asarray([i.masses['gas'] for i in sim.galaxies])
gHI = np.asarray([i.masses['HI'] for i in sim.galaxies])
"""
ii = 0

"""myGal = 7862
gals = np.asarray([i for i in sim.galaxies if i.GroupID==myGal])
"""
for g in gals:
    print 'ii - gmass:',ii,gmass[ii]
    # gas component
    """mstar = np.array([pmass[k] for k in g.glist])
    posstar = np.array([ppos[k] for k in g.glist])
    velstar = np.array([pvel[k] for k in g.glist])
    h1star = np.array([h1[k] for k in g.glist])
    print 'h1star:',np.sum(h1star*mstar),gHI[i]
    """
    mstar = np.array([pmass[k] for k in g.slist])
    posstar = np.array([ppos[k] for k in g.slist])
    velstar = np.array([pvel[k] for k in g.slist])

    posgas = np.array([gaspos[k] for k in g.glist])

    if (len(mstar) > 0):
        min_x = np.min(posstar[:,0])
        max_x = np.max(posstar[:,0])
        min_y = np.min(posstar[:,1])
        max_y = np.max(posstar[:,1])
        min_z = np.min(posstar[:,2])
        max_z = np.max(posstar[:,2])

        binx = np.int((max_x - min_x)/dx)
        biny = np.int((max_y - min_y)/dx)
        binz = np.int((max_z - min_z)/dx)
        print 'bin x,y,z:',binx,biny,binz 

        print 'getting the new center of mass'
        cm = get_cm_radius(posstar,mstar,dx)

        if len(cm) == 4:
            x0,y0,z0 = get_cm(mstar,posstar)
            mask = ( ( (posstar[:,0] - x0)**2 + (posstar[:,1] - y0)**2 + (posstar[:,2] - z0)**2) >= 0.)
            if len(posgas) > 0:
                mask_gas = ( ( (posgas[:,0] - x0)**2 + (posgas[:,1] - y0)**2 + (posgas[:,2] - z0)**2) >= 0.)

        elif len(cm) == 8:
            x0,y0,z0,r0, x2,y2,z2,r2 = cm
            mask = ( ( (posstar[:,0] - x0)**2 + (posstar[:,1] - y0)**2 + (posstar[:,2] - z0)**2) <= r0**2.)
            mask2 = ( ( (posstar[:,0] - x2)**2 + (posstar[:,1] - y2)**2 + (posstar[:,2] - z2)**2) <= r2**2.)
            if len(posgas) > 0:
                mask_gas = (( (posgas[:,0] - x0)**2 + (posgas[:,1] - y0)**2 + (posgas[:,2] - z0)**2) <= r0**2.)
                mask2_gas = (( (posgas[:,0] - x2)**2 + (posgas[:,1] - y2)**2 + (posgas[:,2] - z2)**2) <= r2**2.)

        nb_part.append(len(mstar[mask]))
        mtot.append(np.sum(mstar[mask]))
        slists.append(g.slist[mask])
        new_pos.append(np.array([x0, y0, z0]))
        new_rad.append(r0)
        gal_ids.append(ii)

        if len(posgas) > 0:
            glists.append(g.glist[mask_gas])
        else:
            glists.append(np.array([]))
    
        Lx, Ly, Lz, vrmean0, vzmean0, vtmean0, sigrmean0, sigzmean0, sigtmean0 = get_ang_mtm(mask, posstar, velstar, mstar)

        Lxgal.append(Lx)
        Lygal.append(Ly)
        Lzgal.append(Lz)

        vrmean.append(vrmean0)
        vzmean.append(vzmean0)
        vtmean.append(vtmean0)

        sigrmean.append(sigrmean0)
        sigzmean.append(sigzmean0)
        sigtmean.append(sigtmean0)

        if len(cm) == 8:

            if len(mstar[mask2]) > 256: # do this only for sufficiently massive satellites
                print 'found second galaxy'

                nb_part.append(len(mstar[mask2]))
                mtot.append(np.sum(mstar[mask2]))
                slists.append(g.slist[mask2])
                new_pos.append(np.array([x2, y2, z2]))
                new_rad.append(r2)
                gal_ids.append(ii)
                
                if len(posgas) > 0:
                    glists.append(g.glist[mask2_gas])
                else:
                    glists.append(np.array([]))

                Lx, Ly, Lz, vrmean0, vzmean0, vtmean0, sigrmean0, sigzmean0, sigtmean0 = get_ang_mtm(mask2, posstar, velstar, mstar)

                Lxgal.append(Lx)
                Lygal.append(Ly)
                Lzgal.append(Lz)

                vrmean.append(vrmean0)
                vzmean.append(vzmean0)
                vtmean.append(vtmean0)

                sigrmean.append(sigrmean0)
                sigzmean.append(sigzmean0)
                sigtmean.append(sigtmean0)

    print '\n'
    t += 1
    ii += 1


###output = './Gals_m100n1024_kinematics_'+str(snap)+'_corrcm.asc'
output = './Gals_'+model+'_kinematics_'+str(snap)+'_corrcm_2.asc'
#output = './Gals_m100n1024_kinematics_'+str(snap)+'_stars_r25.asc'
#output = './Gals_m100n1024_kinematics_'+str(snap)+'_stars_r50.asc'
#output = './Gals_m100n1024_kinematics_'+str(snap)+'_stars_r80.asc'
#output = './Gals_m100n1024_kinematics_'+str(snap)+'_gas_r50.asc'
#output = './Gals_m100n1024_kinematics_'+str(snap)+'_gas_r80.asc'

print 'saving into file'

print 'gpos[0]:',gpos[0]
print 'gpos[0,0]:',gpos[:,0][0]

print 'galID=4:',galID[5],Lxgal[5],Lygal[5],Lzgal[5],vtmean[5],sigrmean[5],sigtmean[5],sigzmean[5]

vtmean = np.asarray(vtmean)
sigrmean = np.asarray(sigrmean)
sigzmean = np.asarray(sigzmean)
sigtmean = np.asarray(sigtmean)

print 'shapes:',vtmean.shape,sigrmean.shape,sigzmean.shape,sigtmean.shape

vsig = np.ones(len(gals))
vsig *=-99.
mask = (vtmean != -99.)*(sigrmean != -99.)*(sigzmean != -99.)*(sigtmean != -99.)*(sigrmean != 0.)*(sigtmean != 0.)*(sigzmean != 0.)
vsig[mask] = vtmean[mask]/np.sqrt((sigrmean[mask]**2 + sigzmean[mask]**2 + sigtmean[mask]**2)/3.)

print 'galID=5:',galID[5],Lxgal[5],Lygal[5],Lzgal[5],vtmean[5],sigrmean[5],sigtmean[5],sigzmean[5],vsig[5]

#vsig = np.asarray(vtmean)/np.sqrt((np.asarray(sigrmean)**2 + np.asarray(sigzmean)**2 + np.asarray(sigtmean)**2)/3.)

print 'lengths:',len(galID),len(gmass),len(mtot),len(Lxgal),len(Lygal),len(Lzgal),len(vtmean),len(sigrmean),len(sigzmean),len(sigtmean),len(vsig)

cc = np.array([galID,gmass,mtot,Lxgal,Lygal,Lzgal,vtmean,sigrmean,sigzmean,sigtmean,vsig,nb_part]).T
np.savetxt(output, cc, fmt='%i %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %i', delimiter='\t',header='galID\tMstar\tMstar_r50\tLx\tLy\tLz\tvt\tsigr\tsigz\tsigt\tvsig\tnbPart')

with h5py.File('.Gals_'+model+'_partlists_'+str(snap)+'.h5', 'a') as f:
    f.create_dataset('slists', data=np.array(slists))
    f.create_dataset('glists', data=np.array(glists))
    f.create_dataset('gal_ids', data=np.array(gal_ids))
    f.create_dataset('pos', data=np.array(new_pos))
    f.create_dataset('radius', data=np.array(new_rad))