# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
sys.path.insert(0, '../')
import mscl.plotting

# ----------------------
FIG_NO = 'X6'
colors = mscl.plotting.set_plotting_style()

# Load in the data.
data = pd.read_csv('../data/csv/rbs_predictions.csv')
sorted_vals = data.sort_values('tir_au')
sorted_vals['rel_expression'] = sorted_vals['tir_au'] / \
    sorted_vals[sorted_vals['id'] == 'SD0']['tir_au'].values

sorted_vals
fig, ax = plt.subplots(1, 2, figsize=(6, 3))
ax[0].set_axis_off()
ax[1].yaxis.grid(False)
ax[1].yaxis.set_ticks([0, 1, 2, 3, 4])
ax[1].yaxis.set_ticklabels(['SD6', 'SD4', 'SD2', 'SD1', 'SD0'], fontsize=8)
ax[1].xaxis.set_tick_params(labelsize=8)
ax[1].set_xlabel('predicted expression\n relative to wild-type', fontsize=8)
ax[1].set_xlim([0, 1.1])
_ = ax[1].plot(sorted_vals.rel_expression, [0, 1, 2, 3, 4], 'o',
               color=colors['red'], ms=5)
for i in range(len(sorted_vals)):
    _ = ax[1].hlines(i, 0, sorted_vals.iloc[i]['rel_expression'],
                     color=colors['red'], lw=2)
plt.tight_layout()
plt.savefig('figS{}.svg'.format(FIG_NO))
