import numpy as np
from astropy.io import ascii
import matplotlib.pyplot as plt 

def tukey_biweight(x, c=5.0, epsilon=1.e-11):
    median = np.nanpercentile(x, '50', axis=0)
    mad = np.nanmedian(np.abs(x - median), axis=0)
    u = (x - median)/(c*mad + epsilon)
    weights = np.zeros(u.shape)
    mask = np.abs(u) < 1.
    weights[mask] = (1. - u[mask]**2)**2
    tukey = np.nansum(x*weights, axis=0) / np.sum(weights, axis=0)
    n = x.shape[0]
    num = np.sqrt(np.nansum(((x - tukey)**2. * (1. - u**2.)**4.)*mask, axis=0))
    den = np.abs(np.nansum(((1. - u )*(1 - 5.*(u**2.)))*mask, axis=0))
    scale = np.sqrt(n)*num / den
    t = 1.96 # 97.5% of t distribution with df = max(0.7(n-1), 1)
    return tukey, scale

def plot_belfiore(ax, sample, colors, mass_c=[9.0,  9.5, 10.0,  10.5, 11., 11.5]):
    if sample == 'sf':
        aa = ascii.read('b18_ssfr_gradients_main_sequence.dat')
    if sample == 'gv':
        aa = ascii.read('b18_ssfr_gradients_GV.dat')
    x=aa['R']
    nmass = len(mass_c)
    for nn in range(0, nmass-1):
        y=np.array(aa['{:.2f}'.format(mass_c[nn])+'-'+'{:.2f}'.format(mass_c[nn+1])+' sSFR'])
        y1=np.array(aa['{:.2f}'.format(mass_c[nn])+'-'+'{:.2f}'.format(mass_c[nn+1])+' sSFR']+\
            aa['{:.2f}'.format(mass_c[nn])+'-'+'{:.2f}'.format(mass_c[nn+1])+' error sSFR'])
        y2=np.array(aa['{:.2f}'.format(mass_c[nn])+'-'+'{:.2f}'.format(mass_c[nn+1])+' sSFR']-\
            aa['{:.2f}'.format(mass_c[nn])+'-'+'{:.2f}'.format(mass_c[nn+1])+' error sSFR'])

        #ax.plot(x,y , color=colors[nn], lw=2., label ='{:.1f}'.format(mass_c[nn])+'-'+'{:.1f}'.format(mass_c[nn+1]))
        ax.plot(x,y , color=colors[nn], lw=1.5)
        ax.fill_between(x, y1, y2, where=y1 >= y2, facecolor=colors[nn], alpha=0.2, edgecolor='none')




