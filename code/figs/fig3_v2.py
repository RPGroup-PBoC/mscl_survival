# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.insert(0, '../../')
import mscl.plotting
import mscl.stats

# -------------------
colors = mscl.plotting.set_plotting_style()
FIG_NO = 3

# Load in the data and isolate to onlly the shock experiments.
shock_data = pd.read_csv('../../data/csv/compiled_data.csv')
# shock_data = data[data['experiment'] == 'shock'].copy()

# Define the bins for the range of channels  and colors
bins = np.arange(0, 1100, 30)
color_dict = {True: colors['green'], False: colors['purple']}
label_dict = {True: 'survival', False: 'death'}
zorder_dict = {True: 1, False: 2}
alpha_dict = {True: 0.9, False: 0.8}

# Group by survival and set the colors.
grouped = shock_data.groupby(['survival'])


# Set up the figure canvas.
fig, ax = plt.subplots(1, 2, figsize=(6, 3))

# Format the axes.
ax[0].set_xlabel('effective channel number', fontsize=8)
ax[1].set_xlabel('effective channel number', fontsize=8)
ax[0].set_ylabel('frequency', fontsize=8)
ax[1].set_ylabel('cumulative distribution', fontsize=8)
fig.text(0.01, 0.95, '(A)', fontsize=8)
fig.text(0.51, 0.95, '(B)', fontsize=8)

for a in ax:
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)

for g, d in grouped:
    # Compute the ECDF.
    x, y = mscl.stats.ecdf(d['effective_channels'])
    print(g, np.min(x))
    _ = ax[0].hist(d['effective_channels'], bins=bins, color=color_dict[g], alpha=alpha_dict[g],
                   edgecolor='k', linewidth=0.75, normed=True, label=label_dict[g],
                   zorder=zorder_dict[g])
    _ = ax[1].plot(x, y, '.', ms=2, color=color_dict[g], label='__nolegend__')
    _ = ax[1].plot([], [], '-', marker='.', color=color_dict[g],
                   label=label_dict[g])

# Shade the 100% observed death range.
ax[0].fill_betweenx(np.linspace(0, 0.011, 300), -10,
                    np.min(shock_data[shock_data['survival']
                                      == True]['effective_channels']),
                    color=colors['light_red'],
                    alpha=0.75)

ax[1].fill_betweenx(np.linspace(-0, 1.01, 300), -10,
                    np.min(shock_data[shock_data['survival']
                                      == True]['effective_channels']),
                    color=colors['light_red'],
                    alpha=0.75)

# Add the legend and clean up the figure.
ax[0].legend(fontsize=8)
ax[1].legend(fontsize=8, loc='lower right')

ax[0].set_ylim([0, 0.011])
ax[1].set_ylim([-0.01, 1.01])
ax[1].set_xlim([-10, 850])
ax[0].set_xlim([-10, 850])
plt.tight_layout()
plt.savefig('../../figs/fig{}_dists.pdf'.format(FIG_NO))
plt.savefig('../../figs/fig{}_dists.png'.format(FIG_NO))
