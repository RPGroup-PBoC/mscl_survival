# -*- coding: utf-8 -*-
# %%
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
sys.path.insert(0, '../../')
import mscl.plotting
import mscl.mcmc
import mscl.stats
colors = mscl.plotting.set_plotting_style()


def compute_mean_sem(df):
    n_tot = len(df)
    n_surv = np.sum(df['survival'].astype(int))
    mean_chan = df['effective_channels'].mean()
    prob = n_surv / n_tot
    sem_chan = df['effective_channels'].std() / np.sqrt(len(df))
    prob_err = n_surv * (n_tot - n_surv) / n_tot**3
    out_dict = {'mean_chan': mean_chan, 'p_survival': prob, 'chan_err': sem_chan,
                'prob_err': np.sqrt(prob_err)}
    return pd.Series(out_dict)


# Load the data with the three shock logistic regression.
data = pd.read_csv('../../data/csv/mscl_survival_data_three_shock.csv')

# Load the samples and the stats.
samples = pd.read_csv('../../data/csv/three_shock_complete_mcmc_traces.csv')
stats = pd.read_csv('../../data/csv/three_shock_complete_mcmc_stats.csv')
pooled_samples = pd.read_csv('../../data/csv/pooled_complete_mcmc_traces.csv')
pooled_stats = pd.read_csv('../../data/csv/pooled_complete_mcmc_stats.csv')
pooled_data = pd.read_csv('../../data/csv/pooled_mscl_survival_data.csv')
indiv_samples = pd.read_csv('../../data/csv/complete_mcmc_shock_rate_idx.csv')
indiv_stats = pd.read_csv('../../data/csv/complete_mcmc_shock_rate_idx_stats.csv')

idx = {0: 'slow', 1: 'medium', 2: 'fast'}
stats
# %%
# Instantiate the complicated figure
fig = plt.figure(figsize=(7, 7))
gs = gridspec.GridSpec(22, 30)

# Add panel labels.
fig.text(0.05, 0.9, '(A)', fontsize=10)
fig.text(0.5, 0.9, '(B)', fontsize=10)
fig.text(0.05, 0.48, '(C)', fontsize=10)
fig.text(0.5, 0.48, '(D)', fontsize=10)
# Top panels
surv_ax1 = fig.add_subplot(gs[0, 0:14])
surv_ax2 = fig.add_subplot(gs[0, 16:])
death_ax1 = fig.add_subplot(gs[8, 0:14])
death_ax2 = fig.add_subplot(gs[8, 16:])
surv_ax3 = fig.add_subplot(gs[12, 0:14])
surv_ax4 = fig.add_subplot(gs[12, 16:])
death_ax3 = fig.add_subplot(gs[21, 0:14])
death_ax4 = fig.add_subplot(gs[21, 16:])
points = [surv_ax1, surv_ax2, surv_ax3, surv_ax4,
          death_ax1, death_ax2, death_ax3, death_ax4]

# Survival prob axes.
prob_ax1 = fig.add_subplot(gs[1:8, 0:14])
prob_ax2 = fig.add_subplot(gs[1:8, 16:])
prob_ax3 = fig.add_subplot(gs[13:21, 0:14])
prob_ax4 = fig.add_subplot(gs[13:21, 16:])
curves = [prob_ax1, prob_ax2, prob_ax3, prob_ax4]

# Format the various axes
for ax in points:
    # ax.set_xticks([])
    ax.set_yticks([])
    ax.set_facecolor('none')
    ax.set_xlim([0, 1E3])

for i, ax in enumerate(curves):
    ax.set_xlim([0, 1E3])
    ax.set_ylim([-0.02, 1.02])
    ax.set_xticklabels([])
    if i % 2 == 1:
        ax.set_ylabel('')
        ax.set_yticklabels([])
    else:
        ax.set_ylabel('survival probablity', fontsize=10)
    if i < 2:
        ax.set_xticklabels([])

bottoms = [death_ax1, death_ax2, death_ax3, death_ax4]
for i, ax in enumerate(bottoms):
    ax.set_xlabel('effective channel number', fontsize=10)

# Set the titles
tops = [surv_ax1, surv_ax2, surv_ax3, surv_ax4]
titles = ['slow shock ($< 0.5$ Hz)', 'intermediate shock (0.5 - 1.0 Hz)',
          'fast shock ($ \geq$ 1.0 Hz)', 'all shocks']

for i, ax in enumerate(tops):
    ax.set_title(
        titles[i], backgroundcolor=colors['pale_yellow'], fontsize=10, y=1.1)
    ax.set_xticklabels([])

#  Plot the survival curves and credible regions.
chan_range = np.logspace(0, 3, 500)
for i in range(3):
    beta_0 = stats[stats['parameter'] ==
                   'beta_0__{}'.format(i)]['median'].values[0]
    beta_1 = stats[stats['parameter'] ==
                   'beta_1__{}'.format(i)]['median'].values[0]
    logit = beta_0 + beta_1 * np.log(chan_range)
    prob = (1 + np.exp(-logit))**-1
    _ = curves[i].plot(chan_range, prob, color=colors['red'], lw=1.5)

    # Set the credible regions.
    cred_region = np.zeros((2, len(chan_range)))
    for j, c in enumerate(chan_range):
        logit = samples['beta_0__{}'.format(
            i)] + samples['beta_1__{}'.format(i)] * np.log(c)
        prob = (1 + np.exp(-logit))**-1
        cred_region[:, j] = mscl.mcmc.compute_hpd(prob, 0.95)

    # Fill the shaded cred regions.
    _ = curves[i].fill_between(
        chan_range, cred_region[0, :], cred_region[1, :], color=colors['light_red'], alpha=0.5)

    # Draw the outlines.
    _ = curves[i].plot(chan_range, cred_region[0, :],
                       color=colors['red'], lw=0.5)
    _ = curves[i].plot(chan_range, cred_region[1, :],
                       color=colors['red'], lw=0.5)

# Plot the data.
rates = ['slow', 'medium', 'fast']
for i, r in enumerate(rates):
    _data = data[data['shock_class'] == r]
    surv = _data[_data['survival'] == True]
    death = _data[_data['survival'] == False]
    ys = np.random.normal(0, 0.1, len(surv))
    yd = np.random.normal(0, 0.1, len(death))
    _ = tops[i].plot(surv['effective_channels'], ys, 'k.', ms=1, alpha=0.75)
    _ = bottoms[i].plot(death['effective_channels'],
                        yd, 'k.', ms=1, alpha=0.75)

    # Plot the binned data.
    bin_width = 50
    binned = mscl.stats.density_binning(
        _data, bin_width=bin_width, groupby='shock_class', input_key='effective_channels', min_cells=20)
    grouped = binned.groupby('bin_number').apply(compute_mean_sem)
    _ = curves[i].errorbar(grouped['mean_chan'], grouped['p_survival'], xerr=grouped['chan_err'],
                           yerr=grouped['prob_err'], color='#4b4b4b', lw=1, linestyle='none',
                           marker='.', ms=4,  label='{} channels/bin'.format(bin_width), zorder=100)

    # Plot the strain data.
    grouped = _data.groupby(['rbs']).apply(compute_mean_sem)
    _ = curves[i].errorbar(grouped['mean_chan'], grouped['p_survival'], xerr=grouped['chan_err'], yerr=grouped['prob_err'],
                           color=colors['blue'], lw=1, ms=4, marker='.', label='1 SD mutant / bin', linestyle='none', zorder=100)

# Plot the pooled data.
surv = data[data['survival'] == True]
death = data[data['survival'] == False]
ys = np.random.normal(0, 0.1, len(surv))
yd = np.random.normal(0, 0.1, len(death))
surv_ax4.plot(surv['effective_channels'], ys, 'k.', ms=1, alpha=0.75)
death_ax4.plot(death['effective_channels'], yd, 'k.', ms=1, alpha=0.75)
logit = pooled_stats[pooled_stats['parameter'] == 'beta_0']['median'].values[0] + \
    pooled_stats[pooled_stats['parameter'] ==
                 'beta_1']['median'].values[0] * np.log(chan_range)
prob = (1 + np.exp(-logit))**-1
_ = prob_ax4.plot(chan_range, prob,
                  color=colors['red'], lw=1.5, label='logistic regression')

# Plot the pooled data credible regions.
cred_region = np.zeros((2, len(chan_range)))
for i, c in enumerate(chan_range):
    logit = pooled_samples['beta_0'] + pooled_samples['beta_1'] * np.log(c)
    prob = (1 + np.exp(-logit))**-1
    cred_region[:, i] = mscl.mcmc.compute_hpd(prob, 0.95)
_ = prob_ax4.fill_between(
    chan_range, cred_region[0, :], cred_region[1, :], color=colors['light_red'], alpha=0.5, label='__nolegend__')
_ = prob_ax4.plot(
    chan_range, cred_region[0, :], color=colors['red'], lw=1, label='__nolegend__')
_ = prob_ax4.plot(
    chan_range, cred_region[1, :], color=colors['red'], lw=1, label='__nolegend__')


# Plot the binned pooled data
pooled_data['idx'] = 1
binned = mscl.stats.density_binning(
    pooled_data, bin_width=bin_width, groupby='idx', input_key='effective_channels', min_cells=20)
grouped = binned.groupby('bin_number').apply(compute_mean_sem)
_ = prob_ax4.errorbar(grouped['mean_chan'], grouped['p_survival'], xerr=grouped['chan_err'],
                      yerr=grouped['prob_err'], linestyle='none', lw=1, ms=4, marker='.', color='#4b4b4b', 
                      label='{} channels / bin'.format(bin_width), zorder=100)

# Plot the strain bins.
grouped = pooled_data.groupby(['rbs']).apply(compute_mean_sem)
_ = prob_ax4.errorbar(grouped['mean_chan'], grouped['p_survival'], xerr=grouped['chan_err'], yerr=grouped['prob_err'],
                           color=colors['blue'], lw=1, ms=4, marker='.', label='1 SD mutant / bin', linestyle='none', zorder=100)

_leg = prob_ax4.legend(fontsize=10, handlelength=0.75)

plt.savefig('../../figs/figS4_alternative_binning.pdf', bbox_inches='tight', dpi=300)
plt.savefig('../../figs/figS4_alternative_binning.png', bbox_inches='tight', dpi=300)


# %%
flow_rates = np.sort(data['flow_rate'].unique())


# Set up the complicated figure.
fig, ax = plt.subplots(3, 4, figsize=(8, 5), sharex=True, sharey=True)

# Properly format and label.
for i in range(4):
    ax[-1, i].set_xlabel('effective channel number', fontsize=8)
    if i < 3:
        ax[i, 0].set_ylabel('survival probability', fontsize=8)
ax = ax.ravel()
for a in ax:
    a.tick_params(labelsize=8)
    a.set_xlim([0, 1000])
    a.set_ylim([0, 1])
ax[-1].axis('off')

# Iterate through each shock group and plot the prediction, credible region, and bins.
for i, r in enumerate(flow_rates):
    # Get the most-likely stats.
    beta_0 = indiv_stats[indiv_stats['parameter'] == 'beta_0__{}'.format(i)]['median'].values[0]
    beta_1 = indiv_stats[indiv_stats['parameter'] == 'beta_1__{}'.format(i)]['median'].values[0]

    # Plot the curve.
    logit = beta_0 + beta_1 * np.log(chan_range)
    prob = (1 + np.exp(-logit))**-1

    _ = ax[i].plot(chan_range, prob, color=colors['red'], lw=1)

    # Add the title for the shock rate.
    _ = ax[i].set_title('{} Hz'.format(r), fontsize=8, y=1.04, backgroundcolor=colors['pale_yellow'])

plt.tight_layout()