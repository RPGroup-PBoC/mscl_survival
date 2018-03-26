# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import sys
sys.path.insert(0, '../../../')
import mscl.plotting
import mscl.stats
color = mscl.plotting.set_plotting_style()
FIG_NO = 4
# Load the dataset.
traces = pd.read_csv(
    '../../data/csv/heirarchical_logistic_regression_traces.csv')
data = pd.read_csv('../../data/csv/mscl_survival_data.csv')

# Split by the shock rate.
slow_data = data[(data['experiment'] == 'shock') &
                 (data['shock_class'] == 'slow')].copy()
fast_data = data[(data['experiment'] == 'shock') &
                 (data['shock_class'] == 'fast')].copy()

#%% Compute the survival probability curves given the logistic regression parameters.
chan_range = np.linspace(0, 1000, 500)
samples = ['fast', 'slow']
prob_survival = {}
cred_regions = {}
for s in samples:
    samp = traces[traces['shock_speed'] == s]
    beta_0, beta_1 = np.median(
        samp['beta_0'].values), np.median(samp['beta_1'].values)
    prob_survival[s] = (1 + np.exp(-beta_0 - beta_1 * chan_range))**-1
    _cred = np.zeros((2, len(chan_range)))

    print(s, (np.log(9) - beta_0) / beta_1)

    # Compute the credible regions.
    beta_0, beta_1 = samp['beta_0'].values, samp['beta_1'].values
    print(s, mscl.mcmc.compute_hpd((np.log(9) - beta_0) / beta_1, 0.95))
    for i, c in enumerate(chan_range):
        _prob = (1 + np.exp(-beta_0 - beta_1 * c))**-1
        _cred[:, i] = mscl.mcmc.compute_hpd(_prob, mass_frac=0.95)
    cred_regions[s] = _cred

# %%  Generate figure with appropriate rug plots.
fig = plt.figure(figsize=(6.5, 3.5))
gs = gridspec.GridSpec(3, 2, height_ratios=[0.5, 8, 0.5])
ax = [plt.subplot(gs[i, j]) for i in range(3) for j in range(2)]
ax[0].set_title('slow shock (< 1.0 Hz)', backgroundcolor=color['pale_yellow'],
                fontsize=8)
ax[1].set_title('fast shock (≥ 1.0 Hz)', backgroundcolor=color['pale_yellow'],
                fontsize=8)

# Set the special axes requirements for the rug plots
for i in (0, 1):
    ax[i].set_axis_off()
for a in ax:
    a.tick_params(labelsize=8)
ax[2].set_xticklabels([])
ax[3].set_xticklabels([])

# Plot the survival and death cells on the appropriate rug plots
for j, exp in enumerate([slow_data, fast_data]):
    pos = j % 2
    for i in (0, 1):
        if i == 0:
            grp = True
            a = ax[0 + pos]
        else:
            grp = False
            a = ax[4 + pos]
        _g = exp[exp['survival'] == grp]
        _y = _g['survival'] - np.random.normal(loc=0, scale=0.01, size=len(_g))
        _ = a.plot(_g['effective_channels'], _y, 'k.', ms=1.5, alpha=0.2)

# Properly format the axes labels.
for i in (0, 1, 4, 5):
    ax[i].set_xlim([0, 800])
    ax[i].set_xticks([0, 200, 400, 600, 800])

    if (i == 4) or (i == 5):
        ax[i].set_xticks([0, 200, 400, 600, 800])
        ax[i].set_yticklabels([])
        ax[i].set_facecolor('#FFFFFF')
        ax[i].set_xlabel('effective channel number', fontsize=8)

# Plot the regression curves.
_ = ax[2].plot(chan_range, prob_survival['slow'], color=color['red'],
               label='logistic regression')

_ = ax[3].plot(chan_range, prob_survival['fast'], color=color['blue'])

# Fill in the credible regions.
_ = ax[2].fill_between(chan_range, cred_regions['slow'][0, :], cred_regions['slow'][1, :],
                       color=color['light_red'], alpha=0.5)
_ = ax[2].plot(chan_range, cred_regions['slow'][0, :],
               color=color['red'], lw=0.75, alpha=0.8)
_ = ax[2].plot(chan_range, cred_regions['slow'][1, :],
               color=color['red'], lw=0.75, alpha=0.8)

_ = ax[3].fill_between(chan_range, cred_regions['fast'][0, :], cred_regions['fast'][1, :],
                       color=color['light_blue'], alpha=0.5)
_ = ax[3].plot(chan_range, cred_regions['fast'][0, :],
               color=color['blue'], lw=0.75, alpha=0.6)
_ = ax[3].plot(chan_range, cred_regions['fast'][1, :],
               color=color['blue'], lw=0.75, alpha=0.6)
# Properly set the limits for the regression curves
for i in (2, 3):
    ax[i].set_ylim([0, 1])
    ax[i].set_xlim([0, 1000])
    ax[i].set_ylabel('survival probability', fontsize=8)
plt.subplots_adjust(hspace=0, wspace=0.25)

# Add the appropriate text labels.
ax[0].text(-0.2, 1.5, '(A)', fontsize=8, transform=ax[0].transAxes)
ax[1].text(-0.2, 1.5, '(B)', fontsize=8, transform=ax[1].transAxes)
plt.savefig('../../figs/fig{}.pdf'.format(FIG_NO), bbox_inches='tight')
plt.savefig('../../figs/fig{}.png'.format(FIG_NO), bbox_inches='tight')
