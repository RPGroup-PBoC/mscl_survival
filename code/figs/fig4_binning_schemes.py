# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import sys
sys.path.insert(0, '../../')
import mscl.plotting
import mscl.stats

# -----------------------------------
colors = mscl.plotting.set_plotting_style()
FIG_NO = 4

# Load in the data set and restrict to shock data.
data = pd.read_csv('../../data/csv/mscl_survival_data.csv')
shock_data = data[data['experiment'] == 'shock']


# Set up the figure canvas
fig, ax = plt.subplots(1, 2, figsize=(6, 3.5))

for a in ax:
    a.set_xscale(log)
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
    a.set_ylabel('survival probability', fontsize=8)
    a.set_xlabel('effective channel number', fontsize=8)

ax[0].set_title('binning by strain', fontsize=8,
                backgroundcolor=colors['pale_yellow'], y=1.05)

ax[1].set_title('binning by channel number', fontsize=8,
                backgroundcolor=colors['pale_yellow'], y=1.05)
fig.text(0.0, 0.95, '(A)', fontsize=8)
fig.text(0.5, 0.95, '(B)', fontsize=8)

# Bin by strain and compute the avg..


def compute_mean_sem(df):
    n_tot = len(df)
    n_surv = np.sum(df['survival'].astype(int))
    mean_chan = df['effective_channels'].mean()
    prob = n_surv / n_tot
    sem_chan = df['effective_channels'].std() / np.sqrt(len(df))
    prob_err = n_surv * (n_tot - n_surv) / n_tot**3
    out_dict = {'mean_chan': mean_chan, 'p_survival': prob, 'chan_err': sem_chan,
                'prob_err': np.sqrt(prob_err)}
    return pd.Series(out_dict)


# Compute the statistics and plot
grouped = shock_data.groupby('rbs').apply(compute_mean_sem)
sorted_vals = grouped.sort_values('mean_chan')

sorted_vals
_ = ax[0].errorbar(sorted_vals['mean_chan'], sorted_vals['p_survival'],
                   xerr=sorted_vals['chan_err'], yerr=sorted_vals['prob_err'],
                   marker='.', ms=4, linewidth=1.5, color=colors['red'])


# Define the channel binning widths
color_range = [colors['light_red'], colors['light_blue'],
               colors['green'], colors['purple']]
sorted_vals = shock_data.sort_values('effective_channels')
bin_widths = [10, 50, 100, 200]
for j, w in enumerate(bin_widths):
    bins = np.arange(0, 850, w)
    mean_chan, chan_err, p_survival, prob_err = [], [], [], []
    for i in range(len(bins) - 1):
        slc = sorted_vals.iloc[bins[i]:bins[i + 1] + 1]
        if len(slc) != 0:
            n_surv = np.sum(slc['survival'].astype(int))
            n_tot = len(slc)
            mean_chan.append(slc['effective_channels'].mean())
            chan_err.append(
                slc['effective_channels'].std() / np.sqrt(len(slc)))
            p_survival.append(n_surv / n_tot)
            prob_err.append(n_surv * (n_tot - n_surv) / n_tot**3)

    _ = ax[1].errorbar(mean_chan, p_survival, xerr=chan_err, yerr=np.sqrt(prob_err),
                       color=color_range[j], lw=1.5, marker='.', ms=4, label=w,
                       zorder=j + 1)


leg = ax[1].legend(fontsize=8, title='channels per bin', ncol=2)
plt.setp(leg.get_title(), fontsize=8)
for a in ax:
    a.set_ylim([-0.05, 1])
plt.tight_layout()
plt.savefig('../../figs/fig4_binning.pdf', bbox_inches='tight')
plt.savefig('../../figs/fig4_binning.png', bbox_inches='tight')
