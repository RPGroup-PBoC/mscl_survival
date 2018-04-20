# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import sys
sys.path.insert(0, '../')
import mscl.plotting
import mscl.stats
import mscl.process
colors = mscl.plotting.set_plotting_style()
FIG_NO = 5

data = pd.read_csv('../data/csv/compiled_data.csv')
shock_data = data[data['class'] == 'shock'].copy()
shock_data.loc[shock_data['flow_rate'] < 1, 'shock_speed'] = 'slow'
shock_data.loc[shock_data['flow_rate'] >= 1, 'shock_speed'] = 'fast'
grouped = shock_data.groupby('shock_speed')

# Instantiate and format the figure.
fig, ax = plt.subplots(2, 2, figsize=(6, 5), sharex=True, sharey=True)
ax = ax.ravel()
for a in ax:
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
ax[0].set_ylabel('survival probability', fontsize=8)
ax[2].set_ylabel('survival probability', fontsize=8)
ax[2].set_xlabel('effective channel number', fontsize=8)
ax[3].set_xlabel('effective channel number', fontsize=8)
axes = {'slow': [ax[0], ax[1]], 'fast': [ax[2], ax[3]]}
ax[0].set_title('binning by events', fontsize=8,
                backgroundcolor=colors['pale_yellow'], y=1.05)
ax[1].set_title('binning by channel number', fontsize=8,
                backgroundcolor=colors['pale_yellow'], y=1.05)
fig.text(-0.01, 0.82, 'slow shock ($<$1.0 Hz)', fontsize=8, rotation='vertical',
         backgroundcolor=colors['pale_yellow'])
fig.text(-0.01, 0.385, 'fast shock ($\geq$1.0 Hz)', fontsize=8, rotation='vertical',
         backgroundcolor=colors['pale_yellow'])

# Perform the binning and plot the probability.
sizes = [10, 50, 100, 200]
color_cycle = {'slow': sns.color_palette('Reds_r', n_colors=5),
               'fast': sns.color_palette('Blues_r', n_colors=5)}
for g, d in grouped:
    for i, ev in enumerate(sizes):
        sorted_d = d.sort_values('chan_per_cell')
        mean_chan, prob = [], []
        for j in range(len(d) + ev):
            vals = sorted_d.iloc[j:j + ev]
            if (len(vals) != 0) & (len(vals) >= ev):
                mean_chan.append(vals['chan_per_cell'].mean())
                prob.append(sum(vals['survival'] == True) / len(vals))
        _ = axes[g][0].plot(mean_chan, prob, '-o', label=ev, lw=1, ms=1,
                            color=color_cycle[g][i])

        bins = np.arange(0, np.max(d['chan_per_cell']) + ev, ev)
        mean_chan, prob = [], []
        for k in range(1, len(bins)):
            slc = d[(d['chan_per_cell'] >= bins[k - 1])
                    & (d['chan_per_cell'] < bins[k])]
            if len(slc) != 0:
                mean_chan.append(np.mean(slc['chan_per_cell']))
                prob.append(slc['survival'].sum() / len(slc))

        _ = axes[g][1].plot(mean_chan, prob, '-o', label=ev, lw=1.5, ms=2,
                            color=color_cycle[g][i])

titles = ['events / bin', 'events / bin', 'channels / bin', 'channels / bin']
locs = ['upper right', 'upper left', 'upper right', 'upper left']
for i, a in enumerate(ax):
    _leg = a.legend(fontsize=8, ncol=2, title=titles[i], handlelength=1,
                    columnspacing=0.5, loc=locs[i])
    plt.setp(_leg.get_title(), fontsize=8)


plt.tight_layout()
plt.savefig('figSX{}.pdf'.format(FIG_NO), bbox_inches='tight')
plt.savefig('figSX{}.png'.format(FIG_NO), bbox_inches='tight')
