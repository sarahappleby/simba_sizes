import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib import cm
import sys

def logify(data, error): # for log of x + 1
    log_error = error / (np.log(10.) * (data + 1.))
    log_data = np.log10(data + 1)
    return log_data, log_error

def plot_spindler(ax, labels=False, colors=['b', 'o']):
    spindler_data = {'low_mass': {}, 'high_mass': {}}
    
    spindler_data['low_mass']['r'] = np.array([0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45])
    spindler_data['low_mass']['data'] = np.array([-0.06, -0.005, -0.1, -0.13, -0.195, -0.22, -0.248, -0.28, -0.283, -0.29, -0.3, -0.28, -0.296, -0.266, -0.255])
    spindler_data['low_mass']['lower'] = np.array([-0.184, -0.134, -0.208, -0.234, -0.284, -0.298, -0.319, -0.35, -0.354, -0.354, -0.364, -0.357, -0.365, -0.34, -0.327])
    spindler_data['low_mass']['upper'] = np.array([0.068, 0.121, 0.18, -0.021, -0.096, -0.13, -0.169, -0.201, -0.206, -0.213, -0.22, -0.214, -0.221, -0.191, -0.181])
    
    spindler_data['low_mass']['yerr_upper'] = np.abs(spindler_data['low_mass']['upper'] - spindler_data['low_mass']['data'])
    spindler_data['low_mass']['yerr_lower'] = np.abs(spindler_data['low_mass']['data'] - spindler_data['low_mass']['lower'])

    spindler_data['high_mass']['r'] = np.array([0.05, 0.15, 0.25, 0.35, 0.45, 0.55, 0.65, 0.75, 0.85, 0.95, 1.05, 1.15, 1.25, 1.35, 1.45])
    spindler_data['high_mass']['data'] = np.array([0.335, 0.27, 0.139, 0.015, -0.11, -0.175, -0.2, -0.206, -0.24, -0.193, -0.178, -0.157, -0.147, -0.177, -0.22])
    spindler_data['high_mass']['lower'] = np.array([-0.026, -0.013, -0.094, -0.175, -0.268, -0.309, -0.324, -0.319, -0.335, -0.291, -0.268, -0.245, -0.238, -0.25, -0.284])
    spindler_data['high_mass']['upper'] = np.array([0.695, 0.553, 0.372, 0.21, 0.05, -0.038, -0.089, -0.094, -0.129, -0.09, -0.089, -0.07, -0.056, -0.1, -0.147])
    spindler_data['high_mass']['yerr_upper'] = np.abs(spindler_data['high_mass']['upper'] - spindler_data['high_mass']['data'])
    spindler_data['high_mass']['yerr_lower'] = np.abs(spindler_data['high_mass']['data'] - spindler_data['high_mass']['lower'])

    low_log_data, low_log_error = logify(spindler_data['low_mass']['data'], np.array([spindler_data['low_mass']['yerr_lower'], spindler_data['low_mass']['yerr_upper']]))
    high_log_data, high_log_error = logify(spindler_data['high_mass']['data'], np.array([spindler_data['high_mass']['yerr_lower'], spindler_data['high_mass']['yerr_upper']]))

    if labels:
        ax.errorbar(spindler_data['high_mass']['r'], high_log_data, yerr=[high_log_error[0], high_log_error[1]], capsize=3, c=colors[1], label=labels[1])
        ax.errorbar(spindler_data['low_mass']['r'], low_log_data, yerr=[low_log_error[0], low_log_error[1]], capsize=3, c=colors[0], label=labels[0])
    else:
        ax.errorbar(spindler_data['high_mass']['r'], high_log_data, yerr=[high_log_error[0], high_log_error[1]], capsize=3, c=colors[1])
        ax.errorbar(spindler_data['low_mass']['r'], low_log_data, yerr=[low_log_error[0], low_log_error[1]], capsize=3, c=colors[0])

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rcParams.update({'font.size': 14})

cmap = cm.get_cmap('viridis')
colors = [cmap(0.85), cmap(0.6), cmap(0.35), cmap(0.15)]
colors = [cmap(0.6), cmap(0.3)]

cmap = cm.get_cmap('plasma')
colors_s18 = [cmap(0.6), cmap(0.3)]

plot_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/plots/'
selection = 'sf'

data_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/data/'
plot_dir = '/home/sapple/simba_sizes/profiles/paper/plotting/plots/'
bin_labels = ['10.0-10.6', '>10.6']
plot_labels = [r'$10.0 < \textrm{log} (M_* / M_{\odot}) < 10.6$',
                r'$ \textrm{log} (M_* / M_{\odot}) > 10.6$']
spindler_labels = [r'$\textbf{S18}; 10.0 < \textrm{log} (M_* / M_{\odot}) < 10.6$',
                r'$\textbf{S18}; \textrm{log} (M_* / M_{\odot}) > 10.6$']
xlabel = r'$ R / R_{half}$'

lines = []
for c in colors:
    lines.append(Line2D([0,1],[0,1],linestyle='-', color=c))

line_cen = Line2D([0,1],[0,1],linestyle='--', color='k')
line_sat = Line2D([0,1],[0,1],linestyle=':', color='k')

fig, ax = plt.subplots(2, 2, figsize=(13, 12))
ax = ax.flatten()

leg1 = ax[0].legend(lines, plot_labels, loc=1, fontsize=12)
leg2 = ax[0].legend([line_cen, line_sat],['centrals', 'satellites'], loc=3, fontsize=12)
ax[0].add_artist(leg1)

#leg1 = ax[1].legend(lines, plot_labels, loc=1, fontsize=12)
#leg2 = ax[1].legend([line_cen, line_sat],['centrals', 'satellites'], loc=3, fontsize=12)
#ax[1].add_artist(leg1)

# plot sSFR:
for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+'sf_rot_ssfr_data.h5') as f:
        no_gals = f['cen_no_gals_'+b].value
        tukey = f['cen_tukey_'+b].value
        small = f['cen_small_scale_'+b].value
        large = f['cen_large_scale_'+b].value
        err = f['cen_err_'+b].value
    if i == 0:
        n = len(np.array(tukey))
        dr = 0.2
        factor = dr*n
        bins = np.arange(0., factor, dr)
        rplot = bins+(dr*0.5)

    if i == 1:
        alpha = 0.15
    else:
        alpha = 0.25

    ax[0].plot(rplot, tukey, color=colors[i], linestyle='--')
    ax[0].fill_between(rplot, tukey - err, tukey + err, color=colors[i], alpha=alpha)

    with h5py.File(data_dir+'sf_rot_ssfr_data.h5') as f:
        no_gals = f['sat_no_gals_'+b].value
        tukey = f['sat_tukey_'+b].value
        small = f['sat_small_scale_'+b].value
        large = f['sat_large_scale_'+b].value
        err = f['sat_err_'+b].value
    ax[0].plot(rplot, tukey, color=colors[i], linestyle=':')
    ax[0].fill_between(rplot, tukey - err, tukey + err, color=colors[i], alpha=alpha)

ax[0].set_ylim(-12., -9.5)
ax[0].set_xlim(0, 3.5)
ax[0].set_xlabel(xlabel)
ax[0].set_ylabel(r'$\textrm{log} (\textrm{sSFR} / \textrm{yr}^{-1})$')

for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+'sf_rot_h1_data.h5') as f:
        no_gals = f['cen_no_gals_'+b].value
        tukey = f['cen_tukey_'+b].value
        small = f['cen_small_scale_'+b].value
        large = f['cen_large_scale_'+b].value
        err = f['cen_err_'+b].value
    if i == 0:
        n = len(np.array(tukey))
        dr = 0.2
        factor = dr*n
        bins = np.arange(0., factor, dr)
        rplot = bins+(dr*0.5)

    if i == 1:
        alpha = 0.15
    else:
        alpha = 0.25

    ax[1].plot(rplot, tukey, color=colors[i], linestyle='--')
    ax[1].fill_between(rplot, tukey - err, tukey + err, color=colors[i], alpha=alpha)

    with h5py.File(data_dir+'sf_rot_h1_data.h5') as f:
        no_gals = f['sat_no_gals_'+b].value
        tukey = f['sat_tukey_'+b].value
        small = f['sat_small_scale_'+b].value
        large = f['sat_large_scale_'+b].value
        err = f['sat_err_'+b].value
    ax[1].plot(rplot, tukey, color=colors[i], linestyle=':')
    ax[1].fill_between(rplot, tukey - err, tukey + err, color=colors[i], alpha=alpha)

ax[1].set_xlim(0, 5.)
ax[1].set_ylim(4.5, )
ax[1].set_xlabel(xlabel)
ax[1].set_ylabel(r'$\textrm{log} (\Sigma_{HI} / M_{\odot}\textrm{kpc}^{-2})$')

# plot sSFR:
ax[2].axhline(0, ls='--', c='k')
for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+'sf_rot_sat_cent_ratio_data.h5') as f:
        ssfr_ratio = f['ssfr_ratio_'+b].value
        ssfr_err = np.abs(f['ssfr_err_'+b].value)
    if i == 0:
        n = len(np.array(ssfr_ratio))
        dr = 0.2
        factor = dr*n
        bins = np.arange(0., factor, dr)
        rplot = bins+(dr*0.5)
   
    if i == 0.:
        for j in list(range(1, len(ssfr_err) -1 )):
            if (ssfr_err[j] > 1.):
                ssfr_err[j] = (ssfr_err[j-1] + ssfr_err[j+1]) * 0.5
    if i == 1:
        alpha = 0.15
    else:
        alpha = 0.25

    ax[2].plot(rplot, ssfr_ratio, color=colors[i], linestyle='-')
    ax[2].fill_between(rplot, ssfr_ratio - ssfr_err, ssfr_ratio + ssfr_err, color=colors[i], alpha=alpha)

plot_spindler(ax[2], colors=colors_s18, labels=spindler_labels)
ax[2].legend(loc=1)
ax[2].set_ylim(-1,1)
ax[2].set_xlim(0, 2)
ax[2].set_xlabel(xlabel)
ax[2].set_ylabel(r'$ \textrm{log} (\textrm{sSFR}_{\textrm{sat}}  / \textrm{sSFR}_{\textrm{cen}})$')

# plot H1:
for i, b in enumerate(bin_labels):
    with h5py.File(data_dir+'sf_rot_sat_cent_ratio_data.h5') as f:
        h1_ratio = f['h1_ratio_'+b].value
        h1_err = np.abs(f['h1_err_'+b].value)
    if i == 1:
        alpha = 0.15
    else:
        alpha = 0.25
    ax[3].plot(rplot, h1_ratio, color=colors[i], linestyle='-', label=plot_labels[i])
    ax[3].fill_between(rplot, h1_ratio - h1_err, h1_ratio + h1_err, color=colors[i], alpha=alpha)

ax[3].axhline(0., ls='--', c='k')
ax[3].set_ylim(-0.4, 0.4)
ax[3].set_xlim(0, 2)
ax[3].set_xlabel(xlabel)
ax[3].set_ylabel(r'$ \textrm{log} (\Sigma_{\textrm{HI, sat}} / \Sigma_{\textrm{HI, cen}})$')
#ax[3].legend(loc=1)


plt.savefig(plot_dir+'sat_cent_both.png')
#plt.clf()
