"""
A script to rotate galaxy particle positions and velocities to align with galaxy 
bulk angular momentum, and make images of galaxies face-on and edge-on. Also makes 
a radial profile and fits a Sersic index for each galaxy.

"""
from projection import *
import numpy as np
import caesar
from readgadget import readsnap
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

model = 'm50n512'
snap_no = '078'

data_dir = '/home/rad/data/'+model+'/s50j7k/'
results_dir = '/home/sapple/simba_pretty/quick_projection/'+model+'_'+snap_no+'/galaxies/'
obj = caesar.load(data_dir+'Groups/'+model+'_'+snap_no+'.hdf5')
snap = data_dir + 'snap_'+model+'_'+snap_no+'.hdf5'

h = 0.68 # hubble constant
xmin = -20.   # [kpc]
xmax = 20.    # [kpc]
Npixels = 100 #512
DR = 1. # kpc 
z = obj.simulation.redshift

masses = np.array([i.masses['stellar'] for i in obj.galaxies])
sfr = np.array([i.sfr for i in obj.galaxies])

mask_mass = (masses > 4e10) & (masses < 6e10)
mask_sfr = (sfr > 0.5) & (sfr < 2.5)

#mw_gal_no = np.arange(0, len(masses))[mask_mass*mask_sfr]
mw_gal_no = np.arange(0, len(masses))[mask_mass]

star_positions = readsnap(snap, 'pos', 'star', suppress=1, units=1) / (1.+z)
star_vels = readsnap(snap, 'vel', 'star', suppress=1, units=0)
star_mass = readsnap(snap, 'mass', 'star', suppress=1, units=1) / (h*10**7.)

names = ['edge_on', 'face_on']
vecs = [np.array([0,1,0]), np.array([0,0,1])]

faceon = np.zeros((len(obj.galaxies), 5))
edgeon = np.zeros((len(obj.galaxies), 5))
b_t = np.zeros(len(obj.galaxies))

for i in mw_gal_no:
        print 'galaxy: ' + str(i)
        slist = obj.galaxies[i].slist
        gal_mass = obj.galaxies[i].masses['stellar'].in_units('Msun')
        gal_pos = obj.galaxies[i].pos.in_units('kpc/h')

        x = star_positions[slist][:, 0] - gal_pos[0].value
        y = star_positions[slist][:, 1] - gal_pos[1].value
        z = star_positions[slist][:, 2] - gal_pos[2].value
        vx = star_vels[slist][:, 0]
        vy = star_vels[slist][:, 1]
        vz = star_vels[slist][:, 2]
        mass = star_mass[slist]

        print 'mass:',np.min(mass),np.max(mass)
        print 'x:',np.min(x),np.max(x)

        # ----------------------
        # get projection
        # ----------------------
        r_max = 100.
        r_rot = 500.

        for v in range(len(vecs)):
                try:
                        vec = vecs[v]
                        name = names[v]
                        # 1. compute center of mass for particles within a given radius and shift all particles 
                        posx,posy,posz,vx,vy,vz = recentre_pos_and_vel(x,y,z,vx,vy,vz,mass,r_max)
                        filter_rad = (np.sqrt(posx**2+posy**2+posz**2) < r_rot)
                        # 2. get axis and angle of rotation
                        axis, angle = compute_rotation_to_vec(posx[filter_rad],posy[filter_rad],posz[filter_rad],vx[filter_rad],vy[filter_rad],vz[filter_rad],mass[filter_rad],vec)

                        print 'axis:',axis
                        print 'angle:',angle

                        # 3. rotate positions and velocities
                        posx,posy,posz = rotate(posx,posy,posz,axis,angle)
                        vx,vy,vz = rotate(vx,vy,vz,axis,angle)

                        print 'min-max posx:',np.min(posx),np.max(posx)
                        print 'min-max posy:',np.min(posy),np.max(posy)
                        print 'min-max posz:',np.min(posz),np.max(posz)

                        # -----------------------
                        # Plot:
                        # -----------------------
                        filter_rad=(posx>xmin)*(posx<xmax)*(posy>xmin)*(posy<xmax) ##*(posz>xmin)*(posz<xmax)
                        im,xedges,yedges=np.histogram2d(posx[filter_rad],posy[filter_rad],bins=(Npixels,Npixels),weights=mass[filter_rad])
                        im=im/((xmax-xmin)/float(Npixels))**2 #gives surface density
                        extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]

                        v_min = np.min(np.log10(im[im>0]))
                        v_max = np.max(np.log10(im[im>0]))
                        print 'min-max im:',v_min,v_max

                        plt.imshow(np.log10(im.transpose()+0.0001),extent=extent, interpolation='nearest',cmap='magma',vmin=v_min,vmax=v_max, origin="lower")
                        plt.title('Mass: ' + str.format("{0:.6g}", float(gal_mass)) + ' Msun')
                        plt.savefig(results_dir + 'gal_'+str(i)+ '_'+ name + '.png')
                        plt.clf()

                        mass_new = mass*10**7 #[filter_rad]

                        r = np.sqrt(posx*posx + posy*posy)
                        theta = np.arctan2(posy, posx)    # [-pi,pi]

                        print 'min-max r:',np.min(r),np.max(r)
                        print 'min-max theta:',np.min(theta),np.max(theta)
                        NR = 30#np.int(np.max(r)/DR) # number of radial bins


                        print 'test,NR:',np.int(np.max(r)/DR),NR
                        surface_density = np.zeros(NR)

                        for j in range(0,NR):
                                mask = (r >= j*DR)*(r < (j+1)*DR)
                                surface_density[j] = np.sum(mass_new[mask])

                                if (j==0):
                                        surface_density[j] /= np.pi*DR*DR
                                else:
                                        surface_density[j] /= (np.pi*(DR*DR*(j+1)*(j+1) - DR*DR*j*j))


                        print 'min-max surface_density:',np.min(surface_density),np.max(surface_density)                

                        r_plot = np.arange(0,NR*DR,DR)
                        print 'len:',len(r_plot),len(surface_density)

                        a_guess = np.max(surface_density)
                        r0_guess = np.round(gal_radius, 4)
                        end = np.where(surface_density >= 5e6)[0][-1]
                        r_plot = r_plot[2:end]
                        surface_density = surface_density[2:end]

                        if len(r_plot) > 8:
                                popt, pcov = curve_fit(fit_exp_vac, r_plot, surface_density, p0=(a_guess,1., a_guess, 1.0, r0_guess), maxfev=30000) 
                                print 'popt exp:',popt  # best fit parameters 
                                print '\n'
                                y_exp = fit_exp_vac(r_plot, *popt)  # compute y values using the fitting formula and the best fit parameters
                                
                                plt.semilogy(r_plot,surface_density,'b+',label='data')  
                                plt.semilogy(r_plot,y_exp,'r-',label='fit')
                                plt.xlabel(r"$ R [kpc]$")
                                plt.ylabel(r"$ \Sigma$")
                                plt.title('Mass: ' + str.format("{0:.6g}", float(gal_mass)) + ' Msun')
                                plt.savefig(results_dir + 'gal_'+str(i)+'_'+name+'_exp_vac_fit.png')
                                plt.clf()

                                if name == 'face_on':
                                        faceon[i] = popt
                                elif name == 'edge_on':
                                        edgeon[i] = popt

                                l_bulge = 8*np.pi *popt[2]*(popt[-1]**2)* np.math.factorial(7)
                                l_total = 2*np.pi *popt[0]*(popt[-1]**2) + 8*np.pi *popt[2]*(popt[-1]**2)* np.math.factorial(7)
                                b_t[i] = l_bulge / l_total
        
                        else:
                                pass
                    except (RuntimeError, ValueError) as e:
                        pass

