import numpy as np 
import matplotlib.pyplot as plt
from scipy.ndimage.filters import gaussian_filter

def make_image(posx, posy, weight, filename, rcirc=None, clabel=None, Npixels=100):
	xmin = -20.
	xmax = 20.
	im,xedges,yedges=np.histogram2d(posx,posy,bins=(Npixels,Npixels),weights=weight)
	im=im/((xmax-xmin)/float(Npixels))**2
	sigma = 0.5*Npixels/(xedges.max()-xedges.min())
	im = gaussian_filter(im, sigma)
	extent = [xedges[0], xedges[-1], yedges[0], yedges[-1]]
	v_min = np.min(np.log10(im[im>0]))
	v_max = np.max(np.log10(im[im>0]))

	plt.imshow(np.log10(im.transpose()+0.0001),extent=extent, interpolation='nearest',cmap='magma',
				vmin=v_min,vmax=v_max, origin="lower")
	if rcirc:
		plt.title('Radius: '+ str(rcirc))
	#plt.plot(0., 0., c='w', ms=15, marker='.')
	plt.colorbar(label=clabel)
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
	return thubble - cosmo.age((1./a) -1.).value

# for normalised profiles in physical units:
def make_profile(n, dr, r, quantity, rhalf):

	surface_density = np.zeros(n)
	for j in range(0, n):
		mask = (r >= j*dr)*(r < (j+1)*dr)
		surface_density[j] = np.sum(quantity[mask])
		if (j==0):
			surface_density[j] /= np.pi*(dr*rhalf)**2
		else:
			surface_density[j] /= np.pi*((dr*(j+1)*rhalf)**2 - (dr*j*rhalf)**2) 
	return surface_density

# for physical profiles:
def real_profile(n, dr, r, quantity):
		surface_density = np.zeros(n)
		for j in range(0, n):
				mask = (r >= j*dr)*(r < (j+1)*dr)
				surface_density[j] = np.sum(quantity[mask])
				if (j==0):
						surface_density[j] /= np.pi*(dr**2)
				else:
						surface_density[j] /= np.pi*((dr*(j+1))**2 - (dr*j)**2)
		return surface_density

def npart_profile(n, dr, r):
    prof = np.zeros(n)
    for j in range(0, n):
        mask = (r >= j*dr)*(r < (j+1)*dr)
        prof[i] = len(r[mask])
    return prof

def hi_profile(r, dr, h1_mass, rhalf, h1_limit):
		profile = []
		radius = []
		j = 0
		stop = False
		while not stop:
			mask = (r >= j*dr) * (r < (j+1)*dr)
			profile.append(np.sum(h1_mass*mask))
			radius.append(j*dr + 0.5*dr)
			if j == 0:
				profile[-1] /= np.pi*(dr**2)
			else:
				profile[-1] /= np.pi* ( (dr*(j+1))**2 - (dr*j)**2)
			if (profile[-1] <= h1_limit)  & (j >=int(round(rhalf))):
				stop = True
			else:
				j += 1
		return profile, radius


def plot_profile(r, profile, filename, ylabel, xlabel='R half *', title='', ylim=None):
	plt.plot(r, profile, linestyle='--', marker='.')
	plt.ylabel(ylabel)
	plt.xlabel(xlabel)
	plt.title(title)
	if ylim:
		plt.ylim(ylim[0], ylim[1])
	plt.xlim(0, ) 
	plt.savefig(filename)
	plt.clf()

def plot_h_profile(r, profile, filename, ylabel, xlabel='R half *', title='', ylim=None):
		lower = 4.5
		higher = 8.5
		fig, ax1 = plt.subplots()
		ax1.plot(r, profile, linestyle='--', marker='.')
		ax2 = ax1.twinx()
		if 'h1' in filename:
			ax1.axhline(6., linestyle='--', c='k')
			convert = 1.24e14
			ax2.set_ylim(np.log10(convert*(10**lower)), np.log10(convert*(10**higher)))
			ax2.set_ylabel(r'$ \textrm{log} (N_{HI} / cm^{-2})$')
		elif 'h2' in filename:
			convert = 0.62e14
			ax2.set_ylim(np.log10(convert*(10**lower)), np.log10(convert*(10**higher)))
			ax2.set_ylabel(r'$ \textrm{log} (N_{H_{2}} / cm^{-2})$')
		
		ax1.set_ylabel(ylabel)
		ax1.set_xlabel(xlabel)
		ax1.set_ylim(lower, higher)
		ax1.set_xlim(0, )
		plt.title(title)
		plt.savefig(filename)
		plt.clf()


def bin_data(samples, bins):
	"""
	Bin samples across a range into bins.
	Args:
		samples (ndarray): array of samples
		bins(ndarray): bin values, indicating lower edge of bin
	Returns:
		binned(list): samples sorted into bins
	"""
	width = bins[1] - bins[0]
	binned = []
	for i in range(len(bins)):
		data = []
		for j in range(len(samples)):
			if (samples[j] > bins[i]) & (samples[j] < bins[i]+width):
				data.append(samples[j])
		if not len(data) == 0: binned.append(data)
		else: binned.append([])
	return binned
