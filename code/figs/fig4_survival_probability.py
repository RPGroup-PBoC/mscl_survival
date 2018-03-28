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
print(beta_0, beta_1)


def logit(beta_0, beta_1, chan):
    return (1 + chan**-beta_1 * np.exp(-beta_0))**-1


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


# Compute the prediction from logistic regression.
chan_range = np.linspace(0, 1000, 500)
p_survival = logit(beta_0, beta_1, chan_range)

cred_region = np.zeros((2, len(chan_range)))
for i, c in enumerate(chan_range):
    mode = logit(traces['beta_0'], traces['beta_1'], c)
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
_ = ax1.plot(chan_range, p_survival,
             color=colors['blue'], lw=1, label='no binning')
_ = ax1.fill_between(chan_range, cred_region[0, :], cred_region[1, :],
                     color=colors['light_blue'])
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

# Plot the binned data.
bin_width = 15
binned_data = mscl.stats.density_binning(data, groupby='experiment', bin_width=bin_width,
                                         input_key='effective_channels', min_cells=15)
grouped = binned_data.groupby('bin_number').apply(compute_mean_sem)
_ = ax1.errorbar(grouped['mean_chan'], grouped['p_survival'], xerr=grouped['chan_err'],
                 yerr=grouped['prob_err'], marker='.', color=colors['purple'], linestyle='none',
                 lw=1, label='channel')


# Compute the statistics and plot
grouped = data.groupby('rbs').apply(compute_mean_sem)
sorted_vals = grouped.sort_values('mean_chan')
ax1.errorbar(sorted_vals['mean_chan'], sorted_vals['p_survival'],
             xerr=sorted_vals['chan_err'], yerr=sorted_vals['prob_err'],
             color=colors['red'], linestyle='none', lw=1, marker='.',
             ms=5, zorder=2, label='RBS')
beta_0, beta_1
ax = [ax0, ax1, ax2]
for a in ax:
    # a.set_xscale('log')
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
    a.set_xlim([1, 800])

for a in [ax0, ax2]:
    a.set_frame_on(False)
    _ = a.set_yticks([])

_ = ax0.set_xticks([])
_ = ax1.set_xticklabels(['', '', '', '', ''])
_ = ax1.set_ylim([-0.01, 1.01])
_ = ax1.legend(fontsize=8, title='binning scheme')
plt.setp(_.get_title(), fontsize=8)
plt.subplots_adjust(hspace=0.01)
plt.savefig('../../figs/fig{}_pooled_logistic.pdf'.format(FIG_NO),
            bbox_inches='tight')
plt.savefig('../../figs/fig{}_pooled_logistic.png'.format(FIG_NO),
            bbox_inches='tight')

# %% Generate the pooled data probability curve.
fig, ax = plt.subplots(1, 1, figsize=(4.5, 3.5), height_ratios=[0.5, 0.8, 0.5])


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

chan_range = np.logspace(0, 3, 500)
samples = ['fast', 'slow']
prob_survival = {}
cred_regions = {}
for s in samples:
    samp = traces[traces['shock_speed'] == s]
    beta_0, beta_1 = np.median(
        samp['beta_0'].values), np.median(samp['beta_1'].values)
    prob_survival[s] = logit(beta_0, beta_1, chan_range)

    _cred = np.zeros((2, len(chan_range)))

    print(s, np.exp((np.log(4) - beta_0) / beta_1))
    # Compute the credible regions.
    beta_0, beta_1 = samp['beta_0'].values, samp['beta_1'].values
    print(s, np.exp(mscl.mcmc.compute_hpd((np.log(8) - beta_0) / beta_1, 0.95)))
    for i, c in enumerate(chan_range):
        _prob = logit(beta_0, beta_1, c)
        _cred[:, i] = mscl.mcmc.compute_hpd(_prob, mass_frac=0.95)
    cred_regions[s] = _cred

# %%  Generate figure with appropriate rug plots.
fig = plt.figure(figsize=(6.5, 3.5))
gs = gridspec.GridSpec(3, 2, height_ratios=[0.5, 8, 0.5])
ax = [plt.subplot(gs[i, j]) for i in range(3) for j in range(2)]
ax[0].set_title('slow shock (< 1.0 Hz)', backgroundcolor=colors['pale_yellow'],
                fontsize=8)
ax[1].set_title('fast shock (≥ 1.0 Hz)', backgroundcolor=colors['pale_yellow'],
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
        _y = _g['survival'] - \
            np.random.normal(loc=0, scale=0.01, size=len(_g))
        _ = a.plot(_g['effective_channels'], _y, 'k.', ms=1.5, alpha=0.2)

# Plot data binned by strain
for i, exp in enumerate([slow_data, fast_data]):
    grouped = exp.groupby('rbs').apply(compute_mean_sem)
    _ = ax[i + 2].errorbar(grouped['mean_chan'], grouped['p_survival'], xerr=grouped['chan_err'],
                           yerr=grouped['prob_err'], color=colors['red'], lw=1, linestyle='none',
                           marker='.', ms=4, zorder=100, label='1 RBS/bin')
    binned = mscl.stats.density_binning(
        exp, bin_width=10, groupby='experiment', input_key='effective_channels', min_cells=15)
    grouped = binned.groupby('bin_number').apply(compute_mean_sem)
    _ = ax[i + 2].errorbar(grouped['mean_chan'], grouped['p_survival'], xerr=grouped['chan_err'],
                           yerr=grouped['prob_err'], color='#4b4b4b', lw=1, linestyle='none',
                           marker='.', ms=4, zorder=99, label='10 channels/bin')

# Properly format the axes labels.
for i in (0, 1, 4, 5):

    if (i == 4) or (i == 5):
        # ax[i].set_xticks([0, 200, 400, 600, 800])
        ax[i].set_yticklabels([])
        ax[i].set_facecolor('#FFFFFF')
        ax[i].set_xlabel('effective channel number', fontsize=8)

# Plot the regression curves.
_ = ax[2].plot(chan_range, prob_survival['slow'], color=colors['purple'],
               label='logistic\nregression')

_ = ax[3].plot(chan_range, prob_survival['fast'], color=colors['blue'],
               label='logistic\nregression')

for i in range(2):
    _leg = ax[i + 2].legend(fontsize=8, loc='lower right')
    plt.setp(_leg.get_title(), fontsize=8)
# Fill in the credible regions.
_ = ax[2].fill_between(chan_range, cred_regions['slow'][0, :], cred_regions['slow'][1, :],
                       color=colors['light_purple'], alpha=0.5)
_ = ax[2].plot(chan_range, cred_regions['slow'][0, :],
               color=colors['purple'], lw=0.75, alpha=0.8)
_ = ax[2].plot(chan_range, cred_regions['slow'][1, :],
               color=colors['purple'], lw=0.75, alpha=0.8)
_ = ax[3].fill_between(chan_range, cred_regions['fast'][0, :], cred_regions['fast'][1, :],
                       color=colors['light_blue'], alpha=0.5)
_ = ax[3].plot(chan_range, cred_regions['fast'][0, :],
               color=colors['blue'], lw=0.75, alpha=0.6)
_ = ax[3].plot(chan_range, cred_regions['fast'][1, :],
               color=colors['blue'], lw=0.75, alpha=0.6)
# Properly set the limits for the regression curves
for i in (2, 3):
    ax[i].set_ylim([0, 1])
    ax[i].set_xlim([0, 1000])
    ax[i].set_ylabel('survival probability', fontsize=8)
plt.subplots_adjust(hspace=0, wspace=0.25)

# Add the appropriate text labels.
ax[0].text(-0.2, 1.5, '(A)', fontsize=8, transform=ax[0].transAxes)
ax[1].text(-0.2, 1.5, '(B)', fontsize=8, transform=ax[1].transAxes)
ax[2].set_ylim([-0.01, 1.01])
ax[3].set_ylim([-0.01, 1.01])
for a in ax:
    a.set_xlim([1, 1E3])

plt.savefig('../../figs/fig{}.pdf'.format(FIG_NO), bbox_inches='tight')
plt.savefig('../../figs/fig{}.png'.format(FIG_NO), bbox_inches='tight')
