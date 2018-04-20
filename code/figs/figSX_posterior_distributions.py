# %%
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import pandas as pd
import sys
sys.path.insert(0, '../../')
import mscl.plotting
import mscl.stats
import scipy.stats
colors = mscl.plotting.set_plotting_style()
FIG_NO = 'X1'

# Load in the MCMC data.
data = pd.read_csv('../../data/csv/complete_mcmc_traces.csv')

# Define model contants.
REP = 6

# Set up the figure canvas and format
fig, ax = plt.subplots(1, 2, figsize=(6, 4), sharey=True)
ax[0].set_xlabel('calibration factor [a.u./channel]')
ax[1].set_xlabel('average cell area [µm$^2$]')
ax[0].set_ylabel('$\propto$ probability')
ax[1].set_ylabel('$\propto$ probability')
fig.text(0.0, 1.05, '(A)', fontsize=8)
fig.text(0.55, 1.05, '(B)', fontsize=8)
for a in ax:
    a.set_yticks([])


# Add a fake legend
ax[0].plot([], [], color=colors['blue'], lw=3, label='replicate parameter posterior probability')
ax[0].plot([], [], color=colors['red'], lw=3, label='hyperparameter posterior probability')
ax[0].legend(bbox_to_anchor=(1.05, 1.15), fontsize=8)
# Evaluate the KDE for the low-level parameters
alpha_range = np.linspace(2000, 6000, 500)
area_range = np.linspace(4, 6.5, 500)
for i in range(6):
    # Evaluate the KDE and normalize.
    alpha_kernel = scipy.stats.gaussian_kde(data['alpha__{}'.format(i)])
    area_kernel = scipy.stats.gaussian_kde(data['avg_A__{}'.format(i)])
    alpha_fit = alpha_kernel(alpha_range)
    area_fit = area_kernel(area_range)
    alpha_fit *= np.max(alpha_fit)**-1
    area_fit *= np.max(area_fit)**-1

    # Plot the distributions.
    _ = ax[0].plot(alpha_range, alpha_fit, color=colors['blue'])
    _ = ax[1].plot(area_range, area_fit, color=colors['blue'])
    _ = ax[0].fill_between(alpha_range, alpha_fit, alpha=0.2, color=colors['light_blue'])
    _ = ax[1].fill_between(area_range, area_fit, alpha=0.2, color=colors['light_blue'])


# Plot the hyper parameters
hyper_alpha_kernel = scipy.stats.gaussian_kde(data['hyper_alpha_mu'])
hyper_area_kernel = scipy.stats.gaussian_kde(data['hyper_A_mu'])
hyper_alpha_fit = hyper_alpha_kernel(alpha_range)
hyper_area_fit = hyper_area_kernel(area_range)
hyper_alpha_fit *= np.max(hyper_alpha_fit)**-1
hyper_area_fit *= np.max(hyper_area_fit)**-1
_ = ax[0].plot(alpha_range, hyper_alpha_fit, color=colors['red'], lw=3) 
_ = ax[0].fill_between(alpha_range, hyper_alpha_fit, color=colors['red'], alpha=0.4, zorder=100) 
_ = ax[1].plot(area_range, hyper_area_fit, color=colors['red'], lw=3) 
_ = ax[1].fill_between(area_range, hyper_area_fit, color=colors['red'], alpha=0.4, zorder=100) 

plt.tight_layout()
plt.savefig('../../figs/figS{}.png'.format(FIG_NO), bbox_inches='tight', dpi=300)
plt.savefig('../../figs/figS{}.pdf'.format(FIG_NO), bbox_inches='tight', dpi=300)