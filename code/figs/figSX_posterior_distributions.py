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

#%%
fig = plt.figure(figsize=(6, 4))
gs = gridspec.GridSpec(10, 10)
ax0 = fig.add_subplot(gs[:, 0:6])
ax1 = fig.add_subplot(gs[:5, 6:])
ax2 = fig.add_subplot(gs[5:, 6:])
ax0.set_ylim([2500, 4100])
ax0.set_xlim([4.75, 6])
fig.text(0, 0.95, '(A)', fontsize=8)
fig.text(0.6, 0.95, '(B)', fontsize=8)
fig.text(0.6, 0.5, '(C)', fontsize=8)
# Plot the KDE if the samples.
_ = sns.kdeplot(data['hyper_A_mu'], data['hyper_alpha_mu'], ax=ax0, shade=True,
cmap=plt.cm.Reds)
ax0.set_ylabel('calibration factor [a.u. / MscL channel]', fontsize=8)
ax0.set_xlabel('average cell area [µm$^2$]', fontsize=8)


# Add a fake legend
ax1.plot([], [], color=colors['blue'], lw=3, label='replicate\n parameter')
ax1.plot([], [], color=colors['red'], lw=3, label='hyper-\nparameter')
ax1.legend(fontsize=8, handlelength=1)

# Evaluate the KDE for the low-level parameters
alpha_range = np.linspace(2000, 8000, 500)
area_range = np.linspace(4, 6.5, 500)
for i in range(6):
    # Evaluate the KDE and normalize.
    alpha_kernel = scipy.stats.gaussian_kde(data['alpha__{}'.format(i)])
    area_kernel = scipy.stats.gaussian_kde(data['avg_A__{}'.format(i)])
    alpha_fit = alpha_kernel(alpha_range)
    area_fit = area_kernel(area_range)
    alpha_fit *= np.sum(alpha_fit)**-1
    area_fit *= np.sum(area_fit)**-1

    # Plot the distributions.
    _ = ax1.plot(alpha_range, alpha_fit, color=colors['blue'])
    _ = ax2.plot(area_range, area_fit, color=colors['blue'])
    _ = ax1.fill_between(alpha_range, alpha_fit, alpha=0.2, color=colors['light_blue'])
    _ = ax2.fill_between(area_range, area_fit, alpha=0.2, color=colors['light_blue'])


# Plot the hyper parameters
hyper_alpha_kernel = scipy.stats.gaussian_kde(data['hyper_alpha_mu'])
hyper_area_kernel = scipy.stats.gaussian_kde(data['hyper_A_mu'])
hyper_alpha_fit = hyper_alpha_kernel(alpha_range)
hyper_area_fit = hyper_area_kernel(area_range)
hyper_alpha_fit *= np.sum(hyper_alpha_fit)**-1
hyper_area_fit *= np.sum(hyper_area_fit)**-1
_ = ax1.plot(alpha_range, hyper_alpha_fit, color=colors['red'], lw=3) 
_ = ax1.fill_between(alpha_range, hyper_alpha_fit, color=colors['red'], alpha=0.4, zorder=100) 
_ = ax2.plot(area_range, hyper_area_fit, color=colors['red'], lw=3) 
_ = ax2.fill_between(area_range, hyper_area_fit, color=colors['red'], alpha=0.4, zorder=100) 
ax1.set_xlabel('calibration factor [a.u./channel]', fontsize=8)
ax2.set_xlabel('average cell area [µm$^2$]', fontsize=8)
ax1.set_ylabel('$\propto$ probability', fontsize=8)
ax2.set_ylabel('$\propto$ probability', fontsize=8)
for a in (ax1, ax2):
    a.set_yticks([])

for a in (ax0, ax1, ax2):
    a.tick_params(labelsize=8)
plt.savefig('../../figs/figS{}.png'.format(FIG_NO), bbox_inches='tight', dpi=300)
plt.savefig('../../figs/figS{}.pdf'.format(FIG_NO), bbox_inches='tight', dpi=300)
plt.tight_layout()