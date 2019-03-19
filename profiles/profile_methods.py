import numpy as np 
import matplotlib.pyplot as plt

def make_image(posx, posy, weight, filename, rcirc=None, Npixels=100):
	xmin = -20.
	xmax = 20.
	im,xedges,yedges=np.histogram2d(posx,posy,bins=(Npixels,Npixels),weights=weight)
	im=im/((xmax-xmin)/float(Npixels))**2
	extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
	v_min = np.min(np.log10(im[im>0]))
	v_max = np.max(np.log10(im[im>0]))

	plt.imshow(np.log10(im.transpose()+0.0001),extent=extent, interpolation='nearest',cmap='magma',
				vmin=v_min,vmax=v_max, origin="lower")
	if rcirc:
		plt.title('Radius: '+ str(rcirc))
	plt.plot(0., 0., c='w', ms=15, marker='.')
	plt.colorbar()
	plt.savefig(filename)
	plt.clf()

def center_of_quantity(quantity, weight):
	if len(quantity.shape) == 2:
		weight =  np.transpose(np.array([weight,]*len(quantity[0])))
	return np.sum(quantity*weight, axis=0) / np.sum(weight, axis=0)

def tage(cosmo,thubble,a):
	"""
	Find age of stars in Gyr from expansion time at time of formation
	"""
	return thubble-cosmo.age(1./a-1).value

def make_profile(NR, DR, r, quantity):
	profile = np.zeros(NR)
	# make profile of total mass of stars with ages 50Myr - 100Myr
	for j in range(0, NR):
		mask = (r >= j*DR)*(r < (j+1)*DR)
		profile[j] = np.sum(quantity[mask])
		if (j==0):
			profile[j] /= np.pi*DR*DR
		else:
			profile[j] /= np.pi*(DR*DR*(j+1)*(j+1) - DR*DR*j*j)
	return profile
