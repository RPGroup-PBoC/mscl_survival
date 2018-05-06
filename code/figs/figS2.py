# %%
# -*- coding: utf-8 -*-
import os
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import sys
import glob
sys.path.insert(0, '../../')
import mscl.stats
import mscl.plotting
colors = mscl.plotting.set_plotting_style()
FIG_NO = 2

# -------------
data = pd.read_csv('../../data/csv/mscl_survival_data.csv')
grouped = data.groupby('survival')

# Load the MLG910 data and compute the area
mlg910_files = glob.glob('../processing/*mlg910*/output/*.csv')
mlg910 = pd.concat([pd.read_csv(f, comment='#') for f in mlg910_files], ignore_index=True)
mlg910.loc[mlg910['date'] == 20180410, 'area'] *= 0.5

# Set up the axis
fig, ax = plt.subplots(2, 2, figsize=((6, 6)))
ax = ax.ravel(0)
for i, a in enumerate(ax):
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
    a.set_xlabel('projected cell area [µm$^2$]', fontsize=8)
    if i >= 2: 
        a.set_yscale('log')
        a.set_ylim([5, 1E4])
ax[0].axis('off')

# Labels and titles
ax[1].set_ylabel('cumulative distribution', fontsize=8)
ax[2].set_title('total channel count', fontsize=8, backgroundcolor='#FFEDCE',
                y=1.05)
ax[3].set_title('effective channel count', fontsize=8,
                backgroundcolor='#FFEDCE', y=1.05)
ax[2].set_ylabel('total channels per cell', fontsize=8)
ax[3].set_ylabel('effective channels per cell', fontsize=8)
fig.text(0, 1, '(A)', fontsize=8)
fig.text(0.5, 1, '(B)', fontsize=8)
fig.text(0, 0.5, '(C)', fontsize=8)
fig.text(0.5, 0.5, '(D)', fontsize=8)

# Plot the cumulative distributions for the standard candle and experimental cells
mlg910_x, mlg910_y = mscl.stats.ecdf(mlg910['area'])
data_x, data_y = mscl.stats.ecdf(data['area'])
_ = ax[1].step(mlg910_x, mlg910_y, color='slategray', lw=1.5, label='MLG910')
_ = ax[1].step(data_x, data_y, color=colors['red'], lw=1.5, label='SD mutants')
_ = ax[1].legend(loc='lower right', fontsize=8)
color_dict = {True: colors['green'], False: colors['purple']}
label_dict = {True: 'survival', False: 'death'}
for g, d in grouped:
    _ = ax[2].plot(d['area'], d['scaled_intensity'] * d['area'] / d['calibration_factor'], '.',
                   color=color_dict[g], ms=4, alpha=0.5, label=label_dict[g])
    _ = ax[3].plot(d['area'], d['effective_channels'], '.', color=color_dict[g],
                   ms=4, alpha=0.5, label=label_dict[g])
_ = ax[2].legend(fontsize=8)

plt.tight_layout()
plt.savefig('../../figs/figS{}.pdf'.format(FIG_NO), bbox_inches='tight')
plt.savefig('../figs/figS{}.png'.format(FIG_NO), bbox_inches='tight', dpi=300)
