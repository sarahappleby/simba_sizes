import numpy as np
import h5py


def absmag_to_lum(absmags):
    # assuming absolute AB magnitudes 
    magab0 = 3631e-26 # W m^-2 Hz^-1
    pc = 3.085e16 # m

    flux = np.ones_like(absmags)*10.
    flux  = np.power(flux, absmags*-0.4)

    return flux*magab0*4.*np.pi*(10*pc)**2.


def compute_rfrac(idir0,mass,pos, frac=0.5):
    idir1 = (idir0+1)%3
    idir2 = (idir0+2)%3
    mtot = sum(mass)

    cent = [sum([(pos[i][idir1]) for i in range(len(mass))])/len(mass),sum([(pos[i][idir2]) for i in range(len(mass))])/len(mass)]  # unweighted centroid 
    dpos1 = np.asarray([pos[i][idir1]-cent[0] for i in range(len(mass))])
    dpos2 = np.asarray([pos[i][idir2]-cent[1] for i in range(len(mass))])

    r2 = dpos1*dpos1+dpos2*dpos2  # radii of SF-gas from centroid 
    sortindex = np.argsort(r2)
    r2sort = r2[sortindex]

    msum = np.cumsum(mass[sortindex])

    for i in range(len(mass)):
        if msum[i] > frac*mtot:
            if r2sort[i] > 100*100:
                rh = mygals[igal].radii['stellar_half_mass']
            else: rh = np.sqrt(r2sort[i])
            break
    return rh,cent

# compile magnitudes in desired bands from pyloser file. 
# mtype choices are absmag/appmag/absmag_nodust/appmag_nodust.  also can be spec, iobjs, A_V, or L_FIR (in units of Lsun), in which case bandlist is not needed.
# bandlist is a list of the names of the desired bands, e.g. ['sdss_r','irac_1,'v'].  these must match the band names in the pyloser file.
def get_mags_romeel(loserfile,mtype,bandlist=None,verbose=False):
    hf = h5py.File(loserfile,'r')
    # read in pyloser parameters
    p = [i for i in hf.attrs.items()]
    params = dict([(p[i][0],p[i][1]) for i in range(len(p))])
    if verbose:
        print('Loaded params into dict:',params)
        print('hdf5 file keys: %s' % hf.keys())
    # get desired quantities
    for i in hf.keys():
        if i==mtype: mydata = list(hf[i])
    hf.close()
    print('Retrieved all data')

    # if A_V or L_FIR is requested, these can be directly returned
    if mtype == 'A_V' or mtype == 'L_FIR' or mtype == 'spec' or mtype == 'iobjs': return params,mydata

    # if specific bands are requested, find them in the pyloser file
    bands = params['bands']
    iband1 = iband2 = -1
    magindex = []
    for j in range(len(bandlist)):
        for i in range(len(bands)):
            if bands[i] == bandlist[j]: magindex.append(i)
            #if bands[i].decode('utf-8') == bandlist[j]: magindex.append(i)
    if len(magindex) != len(bandlist):
        print('Bands [',bandlist,'] not all found in pyloser file; available bands=',bands)
        return
    mags = []
    # compile desired band magnitudes
    for j in range(len(bandlist)):
        mags.append(np.asarray([m[magindex[j]] for m in mydata]))
    return params,mags


def get_mags_band(loserfile, mtype, band, only_mtype=False,verbose=False):
    hf = h5py.File(loserfile,'r')
        # read in pyloser parameters
    p = [i for i in hf.attrs.items()]
    params = dict([(p[i][0],p[i][1]) for i in range(len(p))])
    if verbose:
        print('Loaded params into dict:',params)
        print('hdf5 file keys: %s' % hf.keys())
    # get desired quantities
    if only_mtype:
        mags = hf[mtype][:]
    else:
        mags = hf[band+'_'+mtype][:]
    hf.close()
    print('Retrieved all data')
    return params, mags

