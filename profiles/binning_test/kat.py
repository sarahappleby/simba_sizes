import h5py
import sys
import numpy as np
import matplotlib.pyplot as plt

data_file = sys.argv[1]

with h5py.File(data_file, 'r') as f:
	pos = f['pos'][:]
	vel = f['vel'][:]
	mass = f['mass'][:]

x = pos[:, 0]; y = pos[:, 1]; z = pos[:, 2]

fig = plt.figure(1,figsize=(6,5))
fig.subplots_adjust(hspace=0.0, left=0.15, right=0.95, wspace=0.0, bottom=0.1, top=0.93)
ax = fig.add_subplot(111)

ax.set_xlabel(r"$ R [kpc]$")
ax.set_ylabel(r"$ \Sigma $")

# ----------------------------------------
# get polar (r,theta,rho)
# ----------------------------------------
r = np.zeros(len(x))
theta = np.zeros(len(x))
rho_polar = np.zeros(len(x))

r = np.sqrt(x*x + y*y)

print 'r:',r
print 'min-max r:',np.min(r),np.max(r)

##NR = 10 # number of radial bins
maxR = 20
DR = 2.  # step in kpc
NR = np.int(maxR/DR)
rhalf = 10.8 # kpc
rhalf = 1
r = r/rhalf
DR = DR/rhalf

print 'DR,NR:',DR,NR

rho_polar = np.zeros(NR)
surface_density = np.zeros(NR)
# ------------------------------
# 1) no density interpolation:
# ------------------------------
print 'len all:',len(mass)
print 'r min,max:',np.min(r),np.max(r)

for i in range(0,NR):

    print '\n'
    mask = (r >= i*DR)*(r < (i+1)*DR)
    rho_polar[i] = np.sum(mass[mask])   # should be mass here
    surface_density[i] = np.sum(rho_polar[i])
    print len(mass[mask]), np.sum(mass[mask])
    if (i==0):
      surface = np.pi*DR*DR
      print surface_density[i]
      surface_density[i] = surface_density[i]/(np.pi*DR*DR*rhalf*rhalf)
      print 'i:',i,len(mass[mask]),np.sum(mass[mask]),DR 
      print np.pi*DR*DR*rhalf*rhalf
      print surface_density[i]
    else:
      surface = np.pi*(DR*DR*(i+1)*(i+1) - DR*DR*i*i)
      surface_density[i] = surface_density[i]/(np.pi*rhalf*rhalf*(DR*DR*(i+1)*(i+1) - DR*DR*i*i))
      print np.pi*rhalf*rhalf*(DR*DR*(i+1)*(i+1) - DR*DR*i*i)
      print surface_density[i]

    print '\n'
    rho_polar[i] = rho_polar[i]/surface  


r_plot = np.arange(0,NR*DR,DR)

#r_plot = r_plot*phys_size/size  # in kpc, physical size

ax.semilogy(r_plot,surface_density)


plt.show() 


