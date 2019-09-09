import h5py
import numpy as np

data_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/data/'
bin_labels = ['10.0-10.5', '10.5-11.0', '>11.0', '>10.5']
select = 'sf'
gals = 'sat'
fields = ['ssfr', 'sfr', 'fmol', 'h1', 'h2', 'sfe', 'fh2']

for field in fields:
    for i, b in enumerate(bin_labels):
        with h5py.File(data_dir+select+'_rot_'+field+'_data.h5', 'r') as f:
            cv_err = f[gals+'_cv_err_'+b].value
            small = f[gals+'_small_scale_'+b].value
        err = np.sqrt(cv_err**2 + small**2)
        with h5py.File(data_dir+select+'_rot_'+field+'_data.h5', 'a') as f:
            f.create_dataset(gals+'_err_'+b, data=np.array(err))
