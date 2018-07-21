# -*- coding: utf-8 -*-
# %%
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
sys.path.insert(0, '../../')
import mscl.plotting
colors = mscl.plotting.set_plotting_style()
data = pd.read_csv('../../data/csv/poolman_fig5b.csv')

# Define the error reported in van den Berg et al 2016 SI
y_err = 0.3

# Instantiate the figure.
fig, ax = plt.subplots(1, 1, figsize=(4, 3))

# Add labels and format axes
ax.set_xlabel('MscL channels per cell', fontsize=8)
ax.set_ylabel('percent survival', fontsize=8)
ax.set_title('van den Berg et al. 2016 | Figure 5b', fontsize=8,
    backgroundcolor=colors['pale_yellow'], y=1.03)
ax.xaxis.set_tick_params(labelsize=8)
ax.yaxis.set_tick_params(labelsize=8)

# Add data and errorbars
_ = ax.errorbar(data['x'], data['y'], y_err * data['y'], fmt='o', 
                color=colors['red'], ms=3, lw=1)
_ = ax.hlines(data['y'], data['xmin'], data['xmax'], color=colors['red'], lw=1)
plt.tight_layout()
plt.savefig('../../figs/figRX_vandenBerg_comparison.pdf')
plt.savefig('../../figs/figRX_vandenBerg_comparison.png', dpi=300)