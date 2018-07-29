# -*- coding: utf-8 -*-
# %%
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sys.path.insert(0, '../../')
import mscl.plotting
colors = mscl.plotting.set_plotting_style()

# Load the two data sets
figS6 = pd.read_csv('../../data/csv/van_den_Berg_2016_figS6B.csv')
figS6['error'] = figS6['Max'] - figS6['mean']
fig5B = pd.read_csv('../../data/csv/poolman_fig5b.csv')

fig, ax = plt.subplots(1, 2, figsize=(6, 3))

# Format the axes
rotation = (10, 0)
for i, a in enumerate(ax):
    a.xaxis.set_tick_params(labelsize=8, rotation=rotation[i])
    a.yaxis.set_tick_params(labelsize=8)
ax[0].set_ylabel('percent survival', fontsize=8)
ax[1].set_xlabel('MscL channels per cell', fontsize=8)
ax[1].set_ylabel('percent survival', fontsize=8)
ax[0].set_title('van den Berg et al. 2016 - Fig. S6B', fontsize=8,
                backgroundcolor=colors['pale_yellow'], y=1.04)
ax[1].set_title('van den Berg et al. 2016 - Fig. 5B', fontsize=8,
                backgroundcolor=colors['pale_yellow'], y=1.04)

# Add the appropriate panel labels
ax[0].text(-0.3, 1.09, '(A)', fontsize=8, transform=ax[0].transAxes)
ax[1].text(-0.25, 1.09, '(B)', fontsize=8, transform=ax[1].transAxes)

_ = ax[0].bar(width=0.5, bottom=0, x=figS6['sample'], height=figS6['mean'], yerr=figS6['error'], color=colors['red'])
_ = ax[1].errorbar(fig5B['x'], fig5B['y'], xerr=fig5B['xmax'] - fig5B['x'],
               yerr=fig5B['y'] * 0.3, lw=1, fmt='.', color=colors['red'],
               capsize=2)
plt.tight_layout()
plt.savefig('../../figs/figRX_vandenberg_comparison.pdf', bbox_inches='tight')
plt.savefig('../../figs/figRX_vandenberg_comparison.png', bbox_inches='tight', dpi=300)
