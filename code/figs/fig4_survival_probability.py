# -*- coding: utf-8 -*-
# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import sys
sys.path.insert(0, '../../../')
import mscl.plotting
import mscl.stats
colors = mscl.plotting.set_plotting_style()
FIG_NO = 4

# Load the dataset.
traces = pd.read_csv(
    '../../data/csv/pooled_generic_logistic_regression_log.csv')
data = pd.read_csv('../../data/csv/mscl_survival_data.csv')
data = data[data['experiment'] == 'shock']

# Extract the modes.
ind = np.argmax(traces['logp'])
beta_0, beta_1 = traces.iloc[ind][['beta_0', 'beta_1']].values

# Compute the prediction from logistic regression.
chan_range = np.logspace(0, 4, 500)
p_survival = (1 + np.exp(-beta_0 - beta_1 * np.log(chan_range)))**-1

cred_region = np.zeros((2, len(chan_range)))
for i, c in enumerate(chan_range):
    mode = (1 + np.exp(-traces['beta_0'] - traces['beta_1'] * np.log(c)))**-1
    cred_region[:, i] = mscl.mcmc.compute_hpd(mode, 0.95)

# %%
fig = plt.figure(figsize=(5.5, 4.5))
gs = gridspec.GridSpec(3, 1, height_ratios=[0.05, 0.9, 0.05])
ax0 = fig.add_subplot(gs[0])
ax1 = fig.add_subplot(gs[1])
ax2 = fig.add_subplot(gs[2])
ax2.set_xlabel('effective channel number', fontsize=8)
ax1.set_ylabel('survival probability', fontsize=8)

# Plot the prediction and credible region
_ = ax1.plot(chan_range, p_survival, color=colors['blue'], lw=1)
_ = ax1.fill_between(chan_range, cred_region[0, :], cred_region[1, :],
                     color=colors['light_blue'], label='logistic regression')
_ = ax1.plot(chan_range, cred_region[0, :],
             color=colors['blue'], lw=0.5)
_ = ax1.plot(chan_range, cred_region[1, :],
             color=colors['blue'], lw=0.5)

# Plot the survivals and deaths.
surv = data[data['survival'] == True]
death = data[data['survival'] == False]

ax0.plot(surv['effective_channels'], np.random.normal(0, 0.01, size=len(surv)),
         'k.', alpha=0.2, ms=2)
ax2.plot(death['effective_channels'], np.random.normal(0, 0.01,
                                                       size=len(death)), 'k.',
         alpha=0.2, ms=2)

# # Plot the binned data.
# bin_width = 100
# bins = np.arange(0, 850, bin_width)
# sorted_vals = data.sort_values('effective_channels')
# mean_chan, prob, chan_err, prob_err = [], [], [], []
# for i in range(len(bins) - 1):
#     slc = data[(data['effective_channels'] >= bins[i]) &
#                (data['effective_channels'] < bins[i + 1])]
#     if len(slc) > 0:
#         n_surv = np.sum(slc['survival'].astype(int))
#         n_tot = len(slc)
#         mean_chan.append(np.mean(slc['effective_channels']))
#         chan_err.append(np.std(slc['effective_channels']) / np.sqrt(len(slc)))
#         prob.append(n_surv / n_tot)
#         prob_err.append(np.sqrt(n_surv * (n_tot - n_surv) / n_tot**3))
#
# _ = ax1.errorbar(mean_chan, prob, xerr=chan_err, yerr=prob_err, marker='.',
#                  color=colors['red'], linestyle='none', lw=1)
#

# Bin by strain and compute the avg..


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


# Compute the statistics and plot
grouped = data.groupby('rbs').apply(compute_mean_sem)
sorted_vals = grouped.sort_values('mean_chan')
ax1.errorbar(sorted_vals['mean_chan'], sorted_vals['p_survival'],
             xerr=sorted_vals['chan_err'], yerr=sorted_vals['prob_err'],
             color=colors['red'], linestyle='none', lw=1.5, marker='.',
             ms=5, zorder=2)

ax = [ax0, ax1, ax2]
for a in ax:
    a.set_xscale('log')
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
    a.set_xlim([1, 1E4])

for a in [ax0, ax2]:
    a.set_frame_on(False)
    _ = a.set_yticks([])

_ = ax0.set_xticks([])
_ = ax1.set_xticklabels(['', '', '', '', ''])
_ = ax2.set_xticks([1, 10, 100, 1000, 1E4])
_ = ax1.set_ylim([-0.01, 1.01])
plt.subplots_adjust(hspace=0.01)
plt.savefig('../../figs/fig{}_pooled_logistic.pdf', bbox_inches='tight')
plt.savefig('../../figs/fig{}_pooled_logistic.png', bbox_inches='tight')

# %% Generate the pooled data probability curve.
fig, ax = plt.subplots(1, 1, figsize=(4.5, 3.5),
                       height_ratios=[0.5, 0.8, 0.5])


ax.set_xscale('log')

chan_range = np.logspace(0, 4, 200)
prediction = (1 + np.exp(-(beta_0 + beta_1 * np.log(chan_range))))**-1

bins = np.arange(0, 850, 100)
mean_chan, prob = [], []
sorted = data.sort_values('effective_channels')
for i in range(len(bins) - 1):
    slc = sorted.iloc[bins[i]:bins[i + 1]]
    if len(slc) > 0:
        mean_chan.append(slc['effective_channels'].mean())
        prob.append(np.sum(slc['survival'].astype(int)) / len(slc))


ax.plot(chan_range, prediction)
ax.plot(mean_chan, prob, '.', ms=5)

#%% Compute the survival probability curves given the logistic regression parameters.
traces = pd.read_csv(
    '../../data/csv/heirarchical_logistic_regression_traces.csv')

# Split by the shock rate.
slow_data = data[(data['experiment'] == 'shock') &
                 (data['shock_class'] == 'slow')].copy()
fast_data = data[(data['experiment'] == 'shock') &
                 (data['shock_class'] == 'fast')].copy()

chan_range = np.logspace(0, 4, 500)
samples = ['fast', 'slow']
prob_survival = {}
cred_regions = {}
for s in samples:
    samp = traces[traces['shock_speed'] == s]
    beta_0, beta_1 = np.median(
        samp['beta_0'].values), np.median(samp['beta_1'].values)
    prob_survival[s] = (
        1 + np.exp(-beta_0 - beta_1 * np.log(chan_range)))**-1
    _cred = np.zeros((2, len(chan_range)))

    print(s, (np.log(9) - beta_0) / beta_1)

    # Compute the credible regions.
    beta_0, beta_1 = samp['beta_0'].values, samp['beta_1'].values
    print(s, mscl.mcmc.compute_hpd((np.log(9) - beta_0) / beta_1, 0.95))
    for i, c in enumerate(chan_range):
        _prob = (1 + np.exp(-beta_0 - beta_1 * np.log(c)))**-1
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
for a in ax:
    a.set_xscale('log')

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
        _y = _g['survival'] - \
            np.random.normal(loc=0, scale=0.01, size=len(_g))
        _ = a.plot(_g['effective_channels'], _y, 'k.', ms=1.5, alpha=0.2)

# Properly format the axes labels.
for i in (0, 1, 4, 5):

    if (i == 4) or (i == 5):
        # ax[i].set_xticks([0, 200, 400, 600, 800])
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
    # ax[i].set_xlim([0, 1000])
    ax[i].set_ylabel('survival probability', fontsize=8)
plt.subplots_adjust(hspace=0, wspace=0.25)

# Add the appropriate text labels.
ax[0].text(-0.2, 1.5, '(A)', fontsize=8, transform=ax[0].transAxes)
ax[1].text(-0.2, 1.5, '(B)', fontsize=8, transform=ax[1].transAxes)
for a in ax:
    a.set_xlim([1, 1E4])
plt.savefig('../../figs/fig{}.pdf'.format(FIG_NO), bbox_inches='tight')
plt.savefig('../../figs/fig{}.png'.format(FIG_NO), bbox_inches='tight')
