# -*- coding: utf-8 -*-
# %%
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
shock_data = pd.read_csv('../../data/csv/mscl_survival_data.csv')

# Define the bins for the range of channels  and colors
bins = np.arange(0, 1100, 30)
color_dict = {True: colors['green'], False: colors['purple']}
label_dict = {True: 'survival', False: 'death'}
zorder_dict = {True: 1, False: 2}
alpha_dict = {True: 0.9, False: 0.65}

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
    # Compute the ECDF
    x, y = mscl.stats.ecdf(d['effective_channels'])
    xmax = np.sort(d['maximum_channels'])
    xmin = np.sort(d['minimum_channels'])
    _ = ax[0].hist(d['effective_channels'], bins=bins, color=color_dict[g], alpha=alpha_dict[g],
                   edgecolor='k', normed=True, label=label_dict[g],
                   zorder=zorder_dict[g], histtype='stepfilled')
    _ = ax[1].plot(x, y, '.', ms=2, color=color_dict[g], label='__nolegend__')
    _ = ax[1].fill_betweenx(y, xmin, xmax, color=color_dict[g], alpha=0.3, zorder=100)
    _ = ax[1].plot(xmin, y, '-', lw=0.5, color=color_dict[g], label='__nolegend__')
    _ = ax[1].plot(xmax, y, '-', lw=0.5, color=color_dict[g], label='__nolegend__')
    _ = ax[1].plot([], [], '-', marker='.', color=color_dict[g],
                   label=label_dict[g])
sorted_d
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
ax[0].legend(fontsize=10)
ax[1].legend(fontsize=10, loc='lower right')
ax[0].set_ylim([0, 0.011])
ax[1].set_ylim([-0.01, 1.01])
ax[1].set_xlim([-10, 850])
ax[0].set_xlim([-10, 850])
plt.tight_layout()
plt.savefig('../../figs/fig{}.pdf'.format(FIG_NO))
plt.savefig('../../figs/fig{}.png'.format(FIG_NO))


#%% Generate the same plots but separated by shock rate.
grouped = shock_data.groupby(['shock_class', 'survival'])

# Set up the figure canvas.
fig, ax = plt.subplots(1, 2, figsize=(6, 3.5))
fig.text(0, 0.93, '(A)', fontsize=10)
fig.text(0.5, 0.93, '(B)', fontsize=10)
axes = {'slow': ax[0], 'fast': ax[1]}
for g, d in grouped:
    # Compute the ecdf.
    y = np.arange(0, len(d), 1) / len(d)
    x_median = np.sort(d['effective_channels'])
    x_min = np.sort(d['minimum_channels'])
    x_max = np.sort(d['maximum_channels'])

    # Plot the ecdf.
    axes[g[0]].plot(x_median, y, '.', ms=2, color=color_dict[g[1]], label='__nolegend__')
    axes[g[0]].fill_betweenx(y, x_min, x_max, color=color_dict[g[1]], alpha=0.3, label='__nolegend__')
    axes[g[0]].plot(x_max, y, '-', lw=1, color=color_dict[g[1]], label=label_dict[g[1]])
    axes[g[0]].plot(x_min, y, '-', lw=1, color=color_dict[g[1]], label='__nolegend__')
    min_surv = np.min(x_median)
    print(g, min_surv)

    if g[1] == 1:
        axes[g[0]].fill_betweenx(np.linspace(0, 1, 300), 0, min_surv , color=colors['light_red'], alpha=0.5)

ax[0].legend()    
for a in ax:
    a.set_xlim([shock_data['effective_channels'].min(), shock_data['effective_channels'].max()])
    a.set_ylim([0, 1])
    a.set_xlabel('effective channels per cell') 
    a.set_ylabel('cumulative distribution')
ax[0].set_title('slow shock ($<$ 1.0 Hz)', backgroundcolor=colors['pale_yellow'], y=1.04, fontsize=10)
ax[1].set_title('fast shock ($\geq$ 1.0 Hz)', backgroundcolor=colors['pale_yellow'], y=1.04, fontsize=10)
plt.tight_layout()
plt.savefig('../../figs/fig{}_shock_separation.pdf'.format(FIG_NO), bbox_inches='tight')
plt.savefig('../../figs/fig{}_shock_separation.png'.format(FIG_NO), bbox_inches='tight')

shock_data[shock_data['calibration_factor'] > 3000].date.unique()