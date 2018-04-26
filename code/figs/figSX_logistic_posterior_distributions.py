# %%
# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sys.path.insert(0, '../../')
import mscl.stats
import mscl.plotting
import scipy.stats
colors = mscl.plotting.set_plotting_style()
FIG_NO = 'X4'
# Load the data
data = pd.read_csv('../../data/csv/complete_mcmc_traces.csv')
fig, ax = plt.subplots(1, 2, figsize=(6, 4))
fig.text(0, 0.95, '(A)', fontsize=8)
fig.text(0.5, 0.95, '(B)', fontsize=8)
for i, a in enumerate(ax.ravel()):
    a.set_ylabel(r'$\propto$ probability', fontsize=8)
    a.set_yticks([])
    if i < 2:
        a.set_xlabel(r'$\beta_0$', fontsize=8)
        a.set_title(titles[i], fontsize=8,
                    backgroundcolor=colors['pale_yellow'], y=1.08)
    else:
        a.set_xlabel(r'$\beta_1$', fontsize=8)

color_key = ['purple', 'blue']
labels = ['slow shock', 'fast shock']
for i in range(2):
    # Evaluate the KDE of the samples.
    beta0_range = np.linspace(-20, 0, 1000)
    beta1_range = np.linspace(0.5, 3, 1000)
    beta0_kernel = scipy.stats.gaussian_kde(data['beta_0__{}'.format(i)])
    beta0_kde = beta0_kernel(beta0_range)
    beta1_kernel = scipy.stats.gaussian_kde(data['beta_1__{}'.format(i)])
    beta1_kde = beta1_kernel(beta1_range)

    _ = ax[0].plot(beta0_range, beta0_kde, color=colors[color_key[i]],
                   lw=2, zorder=i+1, label=labels[i])
    _ = ax[0].fill_between(beta0_range, beta0_kde,
                           color=colors['light_' + color_key[i]], zorder=i+1, alpha=0.5)
    _ = ax[1].plot(beta1_range, beta1_kde, color=colors[color_key[i]],
                   lw=2, zorder=i+1, label=labels[i])
    _ = ax[1].fill_between(beta1_range, beta1_kde,
                           color=colors['light_' + color_key[i]], zorder=i+1, alpha=0.5)

_ = ax[0].legend(fontsize=8, handlelength=1)
plt.tight_layout()
plt.savefig('../../figs/figS{}_logistic_posterior_distributions.png'.format(FIG_NO),
            bbox_inches='tight')
plt.savefig('../../figs/figS{}_logistic_posterior_distributions.pdf'.format(FIG_NO),
            bbox_inches='tight')
