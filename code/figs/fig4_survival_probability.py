# -*- coding: utf-8 -*-
# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import seaborn as sns
import glob
import sys
sys.path.insert(0, '../../')
import mscl.plotting
import mscl.stats
colors = mscl.plotting.set_plotting_style()
FIG_NO = 4

MAX_EXP = 100
# # Load the dataset.
data = pd.read_csv('../../data/csv/mscl_survival_data.csv')

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

#%% Compute the survival probability curves given the logistic regression parameters.
traces = pd.read_csv(
    '../../data/csv/complete_mcmc_traces.csv')

# Split by the shock rate.
slow_data = data[(data['shock_class'] == 'slow')].copy()
fast_data = data[(data['shock_class'] == 'fast')].copy()

chan_range = np.linspace(0, 1000, 800)
samples = {'fast':1, 'slow':0}
prob_survival = {}
cred_regions = {}
for s in samples.keys():
    beta_0 = traces['beta_0__{}'.format(samples[s])]
    beta_1 = traces['beta_1__{}'.format(samples[s])]
    prob_survival[s] = logit(np.median(beta_0), np.median(beta_1), chan_range)
    _cred = np.zeros((2, len(chan_range)))
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
        _ = a.plot(_g['effective_channels'], _y,  'k.',  ms=1.5, alpha=0.2)
        
# Plot data binned by strain
bin_width=50
for i, exp in enumerate([slow_data, fast_data]):
    grouped = exp.groupby('rbs').apply(compute_mean_sem)
    _ = ax[i + 2].errorbar(grouped['mean_chan'], grouped['p_survival'], xerr=grouped['chan_err'],
                           yerr=grouped['prob_err'], color=colors['red'], lw=1, linestyle='none',
                           marker='.', ms=4, zorder=100, label='1 SD/bin')
    binned = mscl.stats.density_binning(
        exp, bin_width=bin_width, groupby='shock_class', input_key='effective_channels', min_cells=25)
    grouped = binned.groupby('bin_number').apply(compute_mean_sem)
    _ = ax[i + 2].errorbar(grouped['mean_chan'], grouped['p_survival'], xerr=grouped['chan_err'],
                           yerr=grouped['prob_err'], color='#4b4b4b', lw=1, linestyle='none',
                           marker='.', ms=4, zorder=99, label='{} channels/bin'.format(bin_width))

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
    _leg = ax[i + 2].legend(fontsize=8, loc='lower right', handlelength=1)
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
