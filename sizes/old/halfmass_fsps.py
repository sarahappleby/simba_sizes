import matplotlib
matplotlib.use('agg')
import caesar
import sys
import os
import gc
import pylab as plt
import numpy as np
from readgadget import *
import fsps
from astropy.cosmology import FlatLambdaCDM
from cosmocalc import cosmocalc
import multiprocessing as mp

sys.path.append('/home/sapple/tools')
from subdivide import *
import plotmedian as pm

def plot_data(z):
	if z < 0.5: # Zhang+17 SDSS
		ms_data = np.linspace(9,12.0,20)
		alpha = 0.24
		beta = 1.33
		gamma = 10.17
		M0 = 6.49e11/0.68**2
		logRe_spiral = np.log10(gamma * (10**ms_data/M0)**alpha * (1+10**ms_data/M0)**(beta-alpha))
		plt.plot(ms_data,logRe_spiral,':',color='b',lw=4,label='SDSS-LTG (Zhang+17)')
		alpha = 0.17
		beta = 0.58
		gamma = 2.24
		M0 = 2.11e10/0.68**2
		logRe_etg = np.log10(gamma * (10**ms_data/M0)**alpha * (1+10**ms_data/M0)**(beta-alpha))
		plt.plot(ms_data,logRe_etg,':',color='r',lw=4,label='SDSS-ETG (Zhang+17)')
		# Kelvin+11 GAMA
		#logRe_spiral = 0.4+0.3*(ms_data-9)
		#plt.plot(ms_data,logRe_spiral,'-',color='k',lw=4,label='GAMA (Baldry+12)')
	if z >= 0.5 and z < 1.5: # vdWel+14 CANDELS
		logms = [9.25,9.75,10.25,10.75,11.25]
		logRe = [0.37,0.48,0.57,0.67,0.82]
		eRelo = [0.26,0.27,0.24,0.20,0.20]
		eRehi = [0.23,0.21,0.20,0.24,0.14]
		#plt.plot(logms,logRe,'o',color='k',ms=6,label='CANDELS LTG (van der Wel+14)')
		plt.errorbar(np.array(logms)-0.05,logRe,lw=1,yerr=[eRelo,eRehi],fmt='bo',ecolor='b',label='CANDELS LTG')
		logms = [9.25,9.75,10.25,10.75,11.25]
		logRe = [0.23,0.20,0.16,0.38,0.70]
		eRelo = [0.25,0.33,0.26,0.23,0.27]
		eRehi = [0.20,0.34,0.27,0.24,0.23]
		#plt.plot(logms,logRe,'o',color='k',ms=6,label='CANDELS ETG (van der Wel+14)')
		plt.errorbar(np.array(logms)+0.05,logRe,lw=1,yerr=[eRelo,eRehi],fmt='ro',ecolor='r',label='CANDELS ETG')
	if z >= 1.5 and z < 2.5:
		# Alcorn+16 CANDELS+ZFOURGE
		ms_data = np.linspace(9,11.5,5)
		re_data = 0.2*(ms_data-10)+0.4
		plt.plot(ms_data,re_data,'-',color='k',lw=4,label='CANDELS+ZFOURGE (Alcorn+16)')
		# vdWel+14 CANDELS
		logms = [9.75,10.25,10.75,11.25]
		logRe = [0.39,0.44,0.47,0.62]
		eRelo = [0.27,0.32,0.39,0.30]
		eRehi = [0.23,0.21,0.24,0.21]
		plt.plot(logms,logRe,'o',color='k',ms=6,label='CANDELS (van der Wel+14)')
		plt.errorbar(logms,logRe,lw=2,yerr=[eRelo,eRehi],color='k')


def compute_rhalf(idir0,mass,pos):
	idir1 = (idir0+1)%3
	idir2 = (idir0+2)%3
	mtot = sum(mass)
	
	cent = [sum([(mass[i]*pos[i][idir1]) for i in range(len(mass))])/mtot,sum([(mass[i]*pos[i][idir2]) for i in range(len(mass))])/mtot]  # centroid of SFR
	dpos1 = np.asarray([pos[i][idir1]-cent[0] for i in range(len(mass))])
	dpos2 = np.asarray([pos[i][idir2]-cent[1] for i in range(len(mass))])
	r2 = dpos1*dpos1+dpos2*dpos2  # radii of SF-gas from centroid 

	sortindex = np.argsort(r2)
	msum = np.cumsum(mass[sortindex])
	for i in range(len(mass)):
		if msum[i] > 0.5*mtot:
			if r2[i] > 100*100: rh = mygals[igal].radii['stellar_half_mass']
			else: rh = np.sqrt(r2[i])
			break
	return rh,cent

def tage_init(cosmo,redshift):
	thubble = cosmo.age(redshift).value
	lz1 = np.arange(np.log10(1+redshift),np.log10(51),0.0001)
	tlz1 = thubble-cosmo.age(10**lz1-1.).value
	return lz1,tlz1


if __name__ == '__main__':
	if len(sys.argv) < 3:
		print 'usage: MODEL SNAP WIND1 WIND2 ...'
		exit()
	
	# define input file
        MODEL = sys.argv[1]
	SNAP = int(sys.argv[2])
	WIND = [sys.argv[3]]
	#plt.rc('text', usetex=True)
	mlim = 8.7
	mmax = 12.5
	boxsize = int(MODEL[1:3])

	size_comp = True
	plotvar = 'sdss_r'
	#if len(sys.argv)==5: plotvar=sys.argv[4]
	#plotvar = 'sfr'

	plot_dir = '/home/sapple/simba_sizes/plots/'
	data_dir = '/home/sapple/simba_sizes/data/'

	if plotvar not in ['mstar','sfr']:
    		try: sdss_bands = fsps.find_filter(plotvar)
    		except: sys.exit('Filter %s not found'%plotvar)
    		print('Doing rhalf in band %s, generating FSPS stellar pop...'%plotvar)
    		spop = fsps.StellarPopulation(zcontinuous=1, sfh=0, logzsol=0.0, dust_type=2, dust2=0.2)
    		print spop.ssp_ages

	# load in input file
	fig,ax = plt.subplots()
	for iwind in range(0,len(WIND)):
		snapfile = '/home/rad/data/%s/%s/snap_%s_%03d.hdf5' % (MODEL,WIND[iwind],MODEL,SNAP)
		infile = '/home/rad/data/%s/%s/Groups/%s_%03d.hdf5' % (MODEL,WIND[iwind],MODEL,SNAP)
		sim = caesar.load(infile,LoadHalo=False)
		redshift = sim.simulation.redshift
		h = sim.simulation.hubble_constant
		
		mygals = sim.galaxies

		# read in galaxy info
		ids = np.asarray([i.GroupID for i in mygals])
		central = np.asarray([i.central for i in mygals])
		ms = np.asarray([i.masses['stellar'] for i in mygals])
		#mbh = np.asarray([i.masses['bh'] for i in mygals])
		sfr = np.asarray([i.sfr for i in mygals])
		gal_pos = np.asarray([i.pos.in_units('kpc/h') for i in mygals])/(1+redshift)/h # pkpc
		
		ssfr = 1.e9*sfr/ms
		ssfr = np.log10(ssfr+10**(-2.5+0.3*redshift))
		ssfrlim = min(ssfr)+0.2
	 
		# read in particle info
		sm = readsnap(snapfile,'mass','star',units=1,suppress=1)/h  # Mo
		sp = readsnap(snapfile,'pos','star',units=1,suppress=1)/(1+redshift)/h # pkpc
		gm = readsnap(snapfile,'mass','gas',units=1,suppress=1)/h # Mo
		gp = readsnap(snapfile,'pos','gas',units=1,suppress=1)/(1+redshift)/h # pkpc
		gsfr = readsnap(snapfile,'sfr','gas',units=1,suppress=1) # Mo/yr

		if not os.path.isfile(data_dir+'MODEL'+'_'+WIND[0]+'_'+str(SNAP)+'_data.h5'):

			# compute rhalf
			rhalf_l = np.zeros([3,len(mygals)])
			rhalf_m = np.zeros([3,len(mygals)])
			if plotvar=='mstar':
				for igal in range(len(mygals)):
					mass = np.array([sm[k] for k in mygals[igal].slist])
					pos = np.array([sp[k] for k in mygals[igal].slist])
					for idir0 in range(3): 
						rhalf_m[idir0][igal],cent = compute_rhalf(idir0,mass,pos)
					if igal%(len(mygals)/20)==0: 
						print('%d logM*= %.3f c= [%.3f,%.3f] rh2d= %.3f %.3f %.3f rh3d= %.3f' % \
						(igal,np.log10(ms[igal]),cent[0],cent[1],rhalf_m[0][igal],rhalf_m[1][igal],rhalf_m[2][igal],mygals[igal].radii['stellar_half_mass']))
			elif plotvar=='sfr':
				for igal in range(len(mygals)):
					mass = np.array([gsfr[k] for k in mygals[igal].glist])  # "mass" (for weighting) is actually now gas SFR
					pos = np.array([gp[k] for k in mygals[igal].glist])
					if sum(mass)==0: 
						rhalf_m[0][igal] = rhalf_m[1][igal] = rhalf_m[2][igal] = 0
						continue
					for idir0 in range(3): 
						rhalf_m[idir0][igal],cent = compute_rhalf(idir0,mass,pos)
					if igal%(len(mygals)/20)==0: 
						print('%d logM*= %.3f c= [%.3f,%.3f] rh2d= %.3f %.3f %.3f rh3d= %.3f' % \
						(igal,np.log10(ms[igal]),cent[0],cent[1],rhalf_m[0][igal],rhalf_m[1][igal],rhalf_m[2][igal],mygals[igal].radii['stellar_half_mass']))
			else:
				# set up cosmology
				cosmo = FlatLambdaCDM(H0=100*sim.simulation.hubble_constant, 
							Om0=sim.simulation.omega_matter, Ob0=sim.simulation.omega_baryon,Tcmb0=2.73)
				lz1,tlz1 = tage_init(cosmo,redshift)
				zform = 1./readsnap(snapfile,'age','star',units=1,suppress=1)-1
				sa = np.interp(np.log10(zform+1),lz1,tlz1)
				smetal = readsnap(snapfile,'Metallicity','star',units=1,suppress=1) #metallicities
				sZall = np.asarray([i[0] for i in smetal])
				for igal in range(len(mygals)):
					mass = np.array([sm[k] for k in mygals[igal].slist])
					age = np.array([sa[k] for k in mygals[igal].slist])
					sZ = np.array([sZall[k] for k in mygals[igal].slist])
					#Zindex = np.array([np.clip((np.abs(spop.zlegend-i)).argmin()+1,1,len(spop.zlegend)+1) for i in sZ])
					lum = np.zeros(len(sZ))
					for i in range(len(sZ)):
						spop.params['logzsol'] = np.log10(sZ[i]+1.e-6)
						mags = spop.get_mags(tage=age[i],bands=[plotvar])
						lum[i] = mass[i] * 0.4 * 10**(-2.5*mags[0])
					print igal,len(sZ),sZ[:3],lum[:3]
					pos = np.array([sp[k] for k in mygals[igal].slist])
					for idir0 in range(3): 
						rhalf_l[idir0][igal],cent = compute_rhalf(idir0,lum,pos)
						rhalf_m[idir0][igal],cent = compute_rhalf(idir0,mass,pos)
					if igal%(len(mygals)/20)==0: 
						print('%d logM*= %.3f c= [%.3f,%.3f] rh2d= %.3f %.3f %.3f rh3d= %.3f' % \
						(igal,np.log10(ms[igal]),cent[0],cent[1],rhalf_l[0][igal],rhalf_l[1][igal],rhalf_l[2][igal],mygals[igal].radii['stellar_half_mass']))
		                        del mass, age, sZ, pos, lum; gc.collect() 
			#print len(ms),len(sigv3d),len(sigmastar)
			rad_l = (rhalf_l[0]+rhalf_l[1]+rhalf_l[2])/3.
			rad_m = (rhalf_m[0]+rhalf_m[1]+rhalf_m[2])/3.	
			
			with h5py.File(data_dir+MODEL+'_'+WIND[0]+'_'+str(SNAP)+'_data.h5', 'w') as hf:
			    	hf.create_dataset('halflight', data=np.array(rad_l))
			    	hf.create_dataset('halfmass', data=np.array(rad_m))

		condition = (central)
		massbin,cvecbin,ebinlo,ebinhi = pm.runningmedian(np.log10(ms[condition]),ssfr[condition])
		cvec = ssfr - np.interp(np.log10(ms),massbin,cvecbin)
		cmap = plt.get_cmap('jet_r')

		# plot half light radius against stellar mass

		print 'Plotting half light radius against stellar mass'

		fig, ax = plt.subplots()
		xvec = np.log10(ms[rad_l>0])
		yvec = np.log10(rad_l[rad_l>0])
		cvec = cvec[rad_l>0]
		pixsize = 3*(xvec/min(xvec)+1)
		im = ax.scatter(xvec, yvec, c=cvec, s=pixsize, lw=0, cmap=cmap, label=MODEL)
		#plt.plot(xvec, yvec, 'o', c='k', ms=2)
		fig.colorbar(im,ax=ax,label=r'$\Delta\log$sSFR')
		if redshift <= 1.5:
			bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(np.log10(ms[ssfr>ssfrlim]),np.log10(rad_l[ssfr>ssfrlim]))
			ax.plot(bin_cent,ymean,'--',lw=2,color='c')
			#ax.errorbar(bin_cent,ymean,yerr=[ysiglo,ysighi],fmt='ro',c='c')
			bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(np.log10(ms[ssfr<ssfrlim]),np.log10(rad_l[ssfr<ssfrlim]))
			ax.plot(bin_cent,ymean,'--',lw=2,color='m')
			#ax.errorbar(bin_cent,ymean,yerr=[ysiglo,ysighi],fmt='ro',c='r')
		else:
			bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(xvec,yvec)
			ax.plot(bin_cent,ymean,'--',lw=2,color='g')
			#im.set_clim(-0.5,0.5)

		plot_data(redshift)
		#plot_data(0,'.')
		plt.annotate('z=%g,%s'%(np.round(redshift,1),WIND[0]), xy=(0.1, 0.9), xycoords='axes fraction',size=16,bbox=dict(boxstyle="round", fc="w"))
		plt.minorticks_on()
		plt.xlim(mlim,mmax)
		plt.ylim(-0.3,1.8-0.3*redshift)
		plt.xlabel(r'$\log\ M_{*}$',fontsize=16)
		plt.ylabel(r'$\log\ R_{half,*}$' ,fontsize=16)
		plt.legend(loc='lower right')
		plt.savefig(plot_dir+'halflight_'+MODEL+'_'+WIND[0]+'_'+str(SNAP)+'.pdf', bbox_inches='tight', format='pdf')
		plt.clf()
	    
	    
		# plot half mass radius against stellar mass

		print 'Plotting half mass radius against stellar mass'

		fig, ax = plt.subplots()
		xvec = np.log10(ms[(rad_l>0) & (rad_m>0)])
		yvec = np.log10(rad_m[(rad_l>0) & (rad_m>0)])
		cvec = cvec[(rad_l>0) & (rad_m>0)]
		gal_pos = gal_pos[(rad_l>0) & (rad_m>0)]
		pixsize = 3*(xvec/min(xvec)+1)
		im = ax.scatter(xvec, yvec, c=cvec, s=pixsize, lw=0, cmap=cmap, label=MODEL)
		#plt.plot(xvec, yvec, 'o', c='k', ms=2)
		fig.colorbar(im,ax=ax,label=r'$\Delta\log$sSFR')
	 
		if redshift <= 1.5:
			bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(np.log10(rad_m[ssfr>ssfrlim]),np.log10(rad_m[ssfr>ssfrlim]))
			ax.plot(bin_cent,ymean,'--',lw=2,color='c')
			#ax.errorbar(bin_cent,ymean,yerr=[ysiglo,ysighi],fmt='ro',c='c')
			bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(np.log10(rad_m[ssfr<ssfrlim]),np.log10(rad_m[ssfr<ssfrlim]))
			ax.plot(bin_cent,ymean,'--',lw=2,color='m')
			#ax.errorbar(bin_cent,ymean,yerr=[ysiglo,ysighi],fmt='ro',c='r')
	    	else:
			bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(xvec,yvec)
			ax.plot(bin_cent,ymean,'--',lw=2,color='g')

		ax.annotate('z=%g,%s'%(np.round(redshift,1),WIND[0]), xy=(0.1, 0.9), xycoords='axes fraction',size=16,bbox=dict(boxstyle="round", fc="w"))
		ax.minorticks_on()
		plt.xlim(-0.3,1.8-0.3*redshift)
		plt.xlabel(r'$\log\ R_{half,*}$' ,fontsize=16)
		plt.ylim(-0.3,1.8-0.3*redshift)
		plt.ylabel(r'$\log\ R_{half,e}$' ,fontsize=16)
		ax.legend(loc='lower right')
		plt.savefig(plot_dir+'halfmass_'+MODEL+'_'+WIND[0]+'_'+str(SNAP)+'.pdf', bbox_inches='tight', format='pdf')
		plt.clf()

		# plot half light radius against half mass radius

		print 'Plotting half light radius against half mass radius'

		fig, ax = plt.subplots()    
		if redshift <= 1.5:
			bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(np.log10(rad_m[ssfr>ssfrlim]),np.log10(rad_l[ssfr>ssfrlim]))
			ax.plot(bin_cent,ymean,'--',lw=2,color='c', label='SFG')
			#ax.errorbar(bin_cent,ymean,yerr=[ysiglo,ysighi],fmt='ro',c='c')
			bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(np.log10(rad_m[ssfr<ssfrlim]),np.log10(rad_l[ssfr<ssfrlim]))
			ax.plot(bin_cent,ymean,'--',lw=2,color='m', label='low sSFR')
			#ax.errorbar(bin_cent,ymean,yerr=[ysiglo,ysighi],fmt='ro',c='r')
	    	else:
			bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(xvec,yvec)
			ax.plot(bin_cent,ymean,'--',lw=2,color='g')
		#im.set_clim(-0.5,0.5)

		r_bins = np.log10(np.logspace(bin_cent[0], 1.7, num=20))
		r_bins, hlr_mean, hlr_var, hlr_per = jackknife(gal_pos, np.log10(rad_l), np.log10(rad_m), boxsize*1e3, r_bins, std=True)
		ax.errorbar(r_bins, hlr_mean, yerr=[np.zeros(len(hlr_var)), hlr_var], marker='.', markersize=5, linestyle='--', c='b', linewidth=1, capsize=2, 
			label='jackknife')
		ax.fill_between(r_bins, hlr_per[2], hlr_per[1], facecolor='k', alpha=0.2, linewidth=1)
		ax.plot(r_bins, r_bins, c='k', linewidth=1)

		ax.annotate('z=%g,%s'%(np.round(redshift,1),WIND[0]), xy=(0.1, 0.9), xycoords='axes fraction',size=16,bbox=dict(boxstyle="round", fc="w"))
		ax.minorticks_on()
		plt.xlim(-0.3,1.8-0.3*redshift)
		plt.xlabel(r'$\log\ R_{half,*}$' ,fontsize=16)
		plt.ylim(-0.3,1.8-0.3*redshift)
		plt.ylabel(r'$\log\ R_{half,e}$' ,fontsize=16)
		ax.legend(loc='lower right')

		plt.savefig(plot_dir+'halflightmass_median_'+MODEL+'_'+WIND[0]+'_'+str(SNAP)+'.pdf', bbox_inches='tight', format='pdf')
		plt.clf()

