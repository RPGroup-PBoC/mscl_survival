# -*- coding: utf-8 -*-
# %%
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
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

idx = {0: 'slow', 1: 'medium', 2: 'fast'}
# Plot the three.
fig, ax = plt.subplots(2, 2)
ax = ax.ravel()
for a in ax:
    a.set_ylim([0, 1])
    a.set_xlim([0, 1000])
chan_range = np.logspace(0, 3, 500)
for i in range(3):
    beta_0 = stats[stats['parameter'] ==
                   'beta_0__{}'.format(i)]['median'].values[0]
    beta_1 = stats[stats['parameter'] ==
                   'beta_1__{}'.format(i)]['median'].values[0]

    logit = beta_0 + beta_1 * np.log(chan_range)
    prob = (1 + np.exp(-logit))**-1
    ax[i].plot(chan_range, prob, '-', color=colors['red'])

    # Compute the credible regions.
    cred_region = np.zeros((2, len(chan_range)))
    for j, c in enumerate(chan_range):
        logit = samples['beta_0__{}'.format(
            i)] + samples['beta_1__{}'.format(i)] * np.log(c)
        prob = (1 + np.exp(-logit))**-1
        cred_region[:, j] = mscl.mcmc.compute_hpd(prob, 0.95)
    ax[i].fill_between(chan_range, cred_region[0, :],
                       cred_region[1, :], alpha=0.5, color=colors['light_red'])

    # Plot the binned data over
    samps = data[data['shock_class'] == idx[i]]

    grouped = samps.groupby(['rbs']).apply(compute_mean_sem)
    grouped_df = pd.DataFrame(grouped).reset_index()
    ax[i].errorbar(grouped_df['mean_chan'], grouped_df['p_survival'], yerr=grouped_df['prob_err'],
                   xerr=grouped_df['chan_err'], fmt='o', lw=1, ms=3, color=colors['red'])
    binned = mscl.stats.density_binning(samps, bin_width=50, groupby='shock_class', input_key='effective_channels', min_cells=20)
    grouped = binned.groupby('bin_number').apply(compute_mean_sem)
    _ = ax[i].errorbar(grouped['mean_chan'], grouped['p_survival'], xerr=grouped['chan_err'],
                           yerr=grouped['prob_err'], color='#4b4b4b', lw=1, linestyle='none',
                           marker='.', ms=4, zorder=99, label='{} channels/bin'.format(bin_width))

plt.tight_layout()
