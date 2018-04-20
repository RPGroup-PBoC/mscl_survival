# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import pandas as pd
import sys
sys.path.insert(0, '../')
import mscl.plotting
import mscl.stats
colors = mscl.plotting.set_plotting_style()
FIG_NO = 'X1'
# Load and group data by shock group.
traces = pd.read_csv('../data/csv/heirarchical_logistic_regression_traces.csv')
alpha_bins = np.linspace(-5, 0, 75)
beta_bins = np.linspace(0, 0.02, 75)
grouped = traces.groupby('shock_speed')

# Instantiate the figure canvas and set the axes.
fig = plt.figure(figsize=(6.5, 3.5))
gs = gridspec.GridSpec(6, 15)
ax0 = fig.add_subplot(gs[0, :6])
ax1 = fig.add_subplot(gs[1:6, 6])
ax2 = fig.add_subplot(gs[1:6, :6])
ax3 = fig.add_subplot(gs[0, 8:14])
ax4 = fig.add_subplot(gs[1:6, 14])
ax5 = fig.add_subplot(gs[1:6, 8:14])

# Format the axes.
dist_axes = [ax0, ax1, ax3, ax4]
for a in dist_axes:
    a.set_frame_on(False)
    a.set_xticks([])
    a.set_yticks([])
    a.set_xlabel('')
    a.set_ylabel('')

for a in [ax2, ax5]:
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
    a.set_xlabel('intercept', fontsize=8)


ax5.set_ylabel('')
ax5.set_yticks([])
ax2.set_ylabel('slope', fontsize=8)

ax0.set_title('slow shock ($<$ 1.0 Hz)', backgroundcolor='#FFEDCE', fontsize=8)
ax3.set_title(r'fast shock ($\geq$ 1.0 Hz)', backgroundcolor='#FFEDCE',
              fontsize=8)
fig.text(0, 0.92, '(A)', fontsize=8)
fig.text(0.50, 0.92, '(B)', fontsize=8)

# Plot the distributions
axes = {'slow': [ax0, ax1, ax2], 'fast': [ax3, ax4, ax5]}
color_dict = {'slow': [colors['red'], colors['light_red']],
              'fast': [colors['blue'], colors['light_blue']]}
kde_colors = {'slow': 'Reds', 'fast': 'Blues'}
for g, d in grouped:
    _ = sns.distplot(d['alpha'], hist=False, color=color_dict[g][0],
                     kde_kws={'shade': True, 'lw': 1}, ax=axes[g][0])
    _ = sns.distplot(d['beta'], hist=False, color=color_dict[g][0],
                     kde_kws={'shade': True, 'lw': 1}, ax=axes[g][1], vertical=True)
    _ = sns.kdeplot(d['alpha'], d['beta'], shade=True, cmap=kde_colors[g],
                    ax=axes[g][2], clip=((-4.5, -0.5), (0.002, 0.014)))

    # Compute the mean and HPD parameters and plot them on the KDE.
    median_alpha = np.median(d['alpha'])
    median_beta = np.median(d['beta'])
    hpd_alpha = mscl.mcmc.compute_hpd(d['alpha'], mass_frac=0.95)
    hpd_beta = mscl.mcmc.compute_hpd(d['beta'], mass_frac=0.95)

    _ = axes[g][2].plot(median_alpha, median_beta,
                        'slategray', marker='o', ms=2)
    _ = axes[g][2].hlines(median_beta, hpd_alpha[0], hpd_alpha[1],
                          color='slategray', lw=1)
    _ = axes[g][2].vlines(median_alpha, hpd_beta[0], hpd_beta[1],
                          color='slategray', lw=1)


ax5.set_ylabel('')
ax5.set_yticks([])
ax2.set_ylabel('slope', fontsize=8)

plt.tight_layout()
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig('figS{}.pdf'.format(FIG_NO))
plt.savefig('figS{}.png'.format(FIG_NO))
