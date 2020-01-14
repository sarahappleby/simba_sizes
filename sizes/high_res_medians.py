import h5py
import caesar
import numpy as np

import sys
sys.path.append('/home/sapple/tools')
import plotmedian as pm


median_file = '/home/sapple/simba_sizes/sizes/data/m25n512_pyloser_medians.h5'

simba_snaps = ['078', '105', '151']
simba_z = [2.0, 1.0, 0.0]

mtype = 'abs'

model = 'm100n1024'
wind = 's50'

mmin = 10.
mmax = 12.5

for i in range(len(simba_z)):
    caesar_data = '/home/sapple/simba_sizes/sizes/data/'+model+'_'+wind+'_'+simba_snaps[i]+'_caesar_data.h5'
    with h5py.File(caesar_data, 'r') as f:
        central = f['central'].value
        ms = f['stellar_mass'].value
        sfr = f['sfr'].value

    if simba_z[i] == 0.:
        rhalf_file = '/home/sapple/simba_sizes/sizes/data/pyloser_sdss_r.h5'
    else:
        rhalf_file = '/home/sapple/simba_sizes/sizes/data/pyloser_v.h5'
    
    with h5py.File(rhalf_file, 'r') as r:
        rhalf = r[mtype+'_'+model+'_'+wind+'_'+simba_snaps[i]].value
    rhalf = np.sum(rhalf, axis=0) / 3
    
    
    ssfr = 1.e9*sfr/ms
    ssfr = np.log10(ssfr)
    ssfrlim = -1.8+0.3*simba_z[i]

    simba_x = np.log10(ms[rhalf>0])
    simba_y = np.log10(rhalf[rhalf>0])

    with h5py.File(median_file, 'a') as f:
        if 'snap_'+simba_snaps[i] not in f:
            f.create_group('snap_'+simba_snaps[i])

    bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(np.log10(ms[ssfr>ssfrlim]),np.log10(rhalf[ssfr>ssfrlim]))
    with h5py.File(median_file, 'a') as f:
        f['snap_'+simba_snaps[i]][mtype+'_sf_median'] = ymean
        f['snap_'+simba_snaps[i]][mtype+'_sf_bin_cent'] = bin_cent
        f['snap_'+simba_snaps[i]][mtype+'_sf_ysiglo'] = ysiglo
        f['snap_'+simba_snaps[i]][mtype+'_sf_ysighi'] = ysighi

    bin_cent,ymean,ysiglo,ysighi = pm.runningmedian(np.log10(ms[ssfr<ssfrlim]),np.log10(rhalf[ssfr<ssfrlim]))
    with h5py.File(median_file, 'a') as f:
        f['snap_'+simba_snaps[i]][mtype+'_q_median'] = ymean
        f['snap_'+simba_snaps[i]][mtype+'_q_bin_cent'] = bin_cent
        f['snap_'+simba_snaps[i]][mtype+'_q_ysiglo'] = ysiglo
        f['snap_'+simba_snaps[i]][mtype+'_q_ysighi'] = ysighi


