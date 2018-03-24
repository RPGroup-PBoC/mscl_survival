# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import sys
sys.path.insert(0, '../')
import mscl.stats
import mscl.plotting
colors = mscl.plotting.set_plotting_style()
FIG_NO = 2

# -------------
data = pd.read_csv('../../data/csv/mscl_survival_data.csv')

# Examine only the shock data.
shock_data = data[data['experiment'] == 'shock']

# Look only at MLG 910 for channel counts.
mlg910 = data[data['rbs'] == 'mlg910']
grouped = shock_data.groupby('survival')

# Set up the axis
fig, ax = plt.subplots(1, 2, figsize=((6, 3)))
for a in ax:
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
    a.set_xlabel('projected cell area [µm$^2$]', fontsize=8)
    a.set_yscale('log')
    a.set_ylim([5, 5E3])


# Labels and titles
ax[0].set_title('total channel count', fontsize=8, backgroundcolor='#FFEDCE',
                y=1.05)
ax[1].set_title('effective channel count', fontsize=8,
                backgroundcolor='#FFEDCE', y=1.05)
ax[0].set_ylabel('total channels per cell', fontsize=8)
ax[1].set_ylabel('effective channels per cell', fontsize=8)
fig.text(0, 0.92, '(A)', fontsize=8)
fig.text(0.5, 0.92, '(B)', fontsize=8)

# Plotting
color_dict = {True: colors['green'], False: colors['purple']}
label_dict = {True: 'survival', False: 'death'}
for g, d in grouped:
    _ = ax[0].plot(d['area'], d['area'] * d['scaled_intensity'] * CAL, '.',
                   color=color_dict[g], ms=4, alpha=0.5, label=label_dict[g])
    _ = ax[1].plot(d['area'], d['effective_channels'], '.', color=color_dict[g],
                   ms=4, alpha=0.5, label=label_dict[g])
_ = ax[0].legend(fontsize=8)

plt.tight_layout()
plt.savefig('../figs/figSX{}.pdf'.format(FIG_NO), bbox_inches='tight')
plt.savefig('../figs/figSX{}.png'.format(FIG_NO), bbox_inches='tight')
