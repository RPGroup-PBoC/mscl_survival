# -*- coding: utf-8 -*-
#%%
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, '../../')
import mscl.plotting
import mscl.stats
colors = mscl.plotting.set_plotting_style()

# Define the data directory.
data_dir = '../../data/csv/'
data = pd.read_csv('{}mscl_survival_data.csv'.format(data_dir))

# Set up the figure canvas.
fig, ax = plt.subplots(2, 2, figsize=(6, 4), sharex=True)
ax[0,0].set_title('slow shock (< 1.0 Hz)', backgroundcolor=colors['pale_yellow'], y=1.05)
ax[0,1].set_title('fast shock ($\geq$ 1.0 Hz)', backgroundcolor=colors['pale_yellow'], y=1.05)
ax[0, 0].set_ylabel('frequency')
ax[0, 1].set_ylabel('frequency')
ax[1, 0].set_ylabel('cumulative distribution')
ax[1, 1].set_ylabel('cumulative distribution')
axes = {'slow': [ax[0, 0], ax[1, 0]], 'fast': [ax[0, 1], ax[1, 1]]}
color_dict = {True:colors['green'], False: colors['purple']}
labels = {True: 'survival', False: 'death'}
grouped = data.groupby(['shock_class', 'survival'])

# Set the bins
bins = np.arange(0, 800, 50)
for g, d in grouped:
    # Generate the cdf.
    x, y = mscl.stats.ecdf(d['effective_channels'])
    x_max, _ = mscl.stats.ecdf(d['maximum_channels'])
    x_min, _ = mscl.stats.ecdf(d['minimum_channels'])

    # Plot the dists.
    _ = axes[g[0]][0].hist(d['effective_channels'], bins=bins, color=color_dict[g[1]],
    alpha=0.5, lw=1, edgecolor='k', normed=True, label=labels[g[1]])
    # normed=True, alpha=0.5, label=g[1])
    _ = axes[g[0]][1].plot(x, y, '.', color=color_dict[g[1]], alpha=0.5,
    markersize=4, label=labels[g[1]])

ax[0, 0].legend(fontsize=8)
ax[1, 0].legend(fontsize=8)
plt.tight_layout()
plt.savefig('../../figs/split_dists.pdf', bbox_inches='tight')


