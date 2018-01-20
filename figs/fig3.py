import numpy as np
import pandas as pd
import seaborn as sns
import sklearn.neighbors
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pymc3 as pm
import mscl.plotting
import mscl.mcmc
color = mscl.plotting.set_plotting_style()

# Load in the data set and prune.
shock_data = pd.read_csv('../data/csv/compiled_shock_data.csv')
cal_data = pd.read_csv('../data/csv/compiled_calibration_data.csv')

# Remove 10sd1 and frag1
shock_data = shock_data[(shock_data['rbs'] != '10sd1') & (shock_data['date'] != 20170525)
                        & (shock_data['date'] != 20170502)]
cal_data = cal_data[(cal_data['rbs'] != 'frag1') &
                    (cal_data['rbs'] != '10sd1')]

# Merge the two data sets together.
shock_data.loc[:, 'class'] = 'shock'
cal_data.loc[:, 'class'] = 'calibration'
data = pd.concat([shock_data, cal_data], ignore_index=True)

# Convert the RBS identifiers to uppercase.
rbs = list(data['rbs'].values)
new_rbs = [r.upper() for r in rbs]
data.loc[:, 'rbs'] = new_rbs

# Ensure the exposure is set correctly.
frac_exp = data['exposure_ms'].max() / data['exposure_ms']
data.loc[:, 'rescaled_intensity'] = (data['intensity'] -
                                     data['mean_bg']) * frac_exp

# Calculate the standard candle.
CHANNEL_NUM = 348
mlg910 = data[data['rbs'] == 'MLG910']
mean_mlg910 = np.mean(mlg910['rescaled_intensity'] * mlg910['area'])
candle = mean_mlg910 / CHANNEL_NUM
mean_area = mlg910['area'].mean()

# Calculate the relative intensity to MLG910
data.loc[:, 'relative_intensity'] = data['rescaled_intensity'] / \
    mlg910['rescaled_intensity'].mean()

# Calculate the per cell estimate from standard area.
data.loc[:, 'chan_per_cell'] = (
    data['rescaled_intensity'] * mean_area) / candle

# Separate the cells by survivors and fatalities.
surv = data[data['survival'] == True]
fatal = data[data['survival'] == False]

# Compute the ECDFs.
surv_x, surv_y = np.sort(surv['chan_per_cell']), np.arange(
    0, len(surv)) / len(surv)
fatal_x, fatal_y = np.sort(fatal['chan_per_cell']), np.arange(
    0, len(fatal)) / len(fatal)

# Compute the Gaussian KDE for both distributions
chan_range = np.linspace(0, 800, 500)
surv_kern = sklearn.neighbors.KernelDensity(kernel='gaussian',
                                            bandwidth=25).fit(surv['chan_per_cell'].values[:, np.newaxis])
surv_kde = surv_kern.score_samples(chan_range[:, np.newaxis])
fatal_kern = sklearn.neighbors.KernelDensity(
    kernel='gaussian', bandwidth=25).fit(fatal['chan_per_cell'].values[:, np.newaxis])
fatal_kde = fatal_kern.score_samples(chan_range[:, np.newaxis])

# Determine the plotting order.
intensity_order = data.groupby(
    ['rbs'])['relative_intensity'].mean().sort_values()[::-1].index
chan_order = data[data['class'] == 'shock'].groupby(
    ['rbs'])['chan_per_cell'].mean().sort_values()[::-1].index
shock_exp = data[data['class'] == 'shock']

#%% Set up the plots.
fig = plt.figure(figsize=(6.5, 5))
gs = gridspec.GridSpec(2, 2)

ax0 = plt.subplot(gs[0, 0])
ax1 = plt.subplot(gs[0, 1])
ax2 = plt.subplot(gs[1, :])

# Plot the intensity data.
_ = sns.boxplot(x='rbs', y='relative_intensity', data=data, order=intensity_order,
                linewidth=0.5, palette='Greens', ax=ax0, fliersize=0)
_ = sns.stripplot(x='rbs', y='relative_intensity', data=data, order=intensity_order,
                  marker='.', color='k', size=1, ax=ax0, jitter=True)

# Add a marker for the candle strain.
ax0.vlines(2, ax0.get_ylim()[0], ax0.get_ylim()[1],
           color=color['pale_yellow'], lw=20, alpha=0.75, zorder=0)

# Plot the shock data boxplots
_ = sns.boxplot(x='rbs', y='chan_per_cell', data=data, order=chan_order,
                linewidth=0.5, palette='Greens_r', ax=ax1, fliersize=0,
                hue='survival')

# Plot the ECDFs and KDE.
surv_colors = sns.color_palette('Greens_r', n_colors=4)
_ = ax2.plot(surv_x, surv_y, '-', color=surv_colors[0], label='survival')
_ = ax2.plot(fatal_x, fatal_y, '-', color=surv_colors[2], label='death')

# Set the labels.
ax0.set_xlabel('RBS identifier', fontsize=8)
ax0.set_ylabel('relative expression', fontsize=8)
ax1.set_xlabel('RBS identifier', fontsize=8)
ax1.set_ylabel('channels per cell', fontsize=8)
ax2.set_xlabel('channels per cell', fontsize=8)
ax2.set_ylabel('ECDF', fontsize=8)

# Set the legend for axis 2
hand, lab = ax1.get_legend_handles_labels()
ax1.legend(hand, ['survival', 'death'], fontsize=8)
ax2.legend(['survival', 'death'], fontsize=8)

# Add panel labels
ax0.text(-0.18, 0.95, '(A)',  transform=ax0.transAxes, fontsize=10)
ax1.text(-0.18, 0.95, '(B)',  transform=ax1.transAxes, fontsize=10)
ax2.text(-0.08, 0.95, '(C)', transform=ax2.transAxes, fontsize=10)
ax2.set_xlim([0, 800])
# Clean up the figure
ax0.set_ylim([0, 2.75])
for t in ax0.get_xticklabels():
    t.set_rotation(45)

for a in [ax0, ax1, ax2]:
    a.tick_params(labelsize=8)
plt.tight_layout()

# Save the figure as PDF and PNG.
plt.savefig('fig3.pdf', bbox_inches='tight', dpi=300)
plt.savefig('fig3.png', bbox_inches='tight', dpi=300)


# %% Figure 4 - Logistic Regression.

# Separate data by slow and fast.
ident = {True: 'fast', False: 'slow'}
group = [ident[i >= 1.0] for i in shock_exp['flow_rate'].values]
shock_exp.loc[:, 'shock_group'] = group

# Perform the logistic regression.
with pm.Model() as slow_model:
    pm.glm.GLM.from_formula('survival ~ chan_per_cell',
                            data=shock_exp[(shock_exp['shock_group'] == 'slow') &
                                           (shock_exp['chan_per_cell'] <= 700)],
                            family=pm.families.Binomial())
    trace = pm.sample(draws=5000, tune=5000, njobs=4)
    slow_trace = mscl.mcmc.trace_to_dataframe(trace, slow_model)
    slow_stats = mscl.mcmc.compute_statistics(slow_trace)


with pm.Model() as fast_model:
    pm.glm.GLM.from_formula('survival ~ chan_per_cell',
                            data=shock_exp[shock_exp['shock_group'] == 'fast'],
                            family=pm.families.Binomial())
    trace = pm.sample(draws=5000, tune=5000, njobs=4)
    fast_trace = mscl.mcmc.trace_to_dataframe(trace, fast_model)
    fast_stats = mscl.mcmc.compute_statistics(fast_trace)

with pm.Model() as two_dim_model:
    pm.glm.GLM.from_formula('survival ~ chan_per_cell + flow_rate',
                            data=shock_exp, family=pm.families.Binomial())
    trace = pm.sample(draws=5000, tune=5000, njobs=4)
    two_dim_trace = mscl.mcmc.trace_to_dataframe(trace, two_dim_model)
    two_dim_stats = mscl.mcmc.compute_statistics(two_dim_trace)

# %%

# Compute the logistic regression curves for the slow and fast groups.
slow_beta_0, slow_beta_1 = slow_stats['mode'].values
fast_beta_0, fast_beta_1 = fast_stats['mode'].values
chan_range = np.linspace(0, 800, 200)
slow_prob = (1 + np.exp(slow_beta_0 + slow_beta_1 * chan_range))**-1
fast_prob = (1 + np.exp(fast_beta_0 + fast_beta_1 * chan_range))**-1

# Compute the credible regions.
slow_cred_region = np.zeros(shape=(2, len(chan_range)))
fast_cred_region = np.zeros(shape=(2, len(chan_range)))
for i, c in enumerate(chan_range):
    # Evaluate the probability.
    _slow_prob = (
        1 + np.exp(slow_trace['Intercept'] + c * slow_trace['chan_per_cell']))**-1
    _fast_prob = (
        1 + np.exp(fast_trace['Intercept'] + c * fast_trace['chan_per_cell']))**-1
    slow_cred_region[:, i] = mscl.mcmc.compute_hpd(_slow_prob, mass_frac=0.95)
    fast_cred_region[:, i] = mscl.mcmc.compute_hpd(_fast_prob, mass_frac=0.95)


# Separate the survivals from the deaths.
slow_exp = shock_exp[(shock_exp['shock_group'] == 'slow')
                     & (shock_exp['chan_per_cell'] <= 700)]
fast_exp = shock_exp[shock_exp['shock_group'] == 'fast']

# Jitter them around 0 and 1.
slow_plot = slow_exp['survival'].values.astype(
    int) - np.random.normal(loc=0, scale=0.01, size=len(slow_exp))
fast_plot = fast_exp['survival'].values.astype(
    int) - np.random.normal(loc=0, scale=0.01, size=len(fast_exp))

# %% Set up the rug plots.
fig, ax = plt.subplots(1, 2, figsize=(6.5, 3))

ax[0].set_xlabel('channels per cell', fontsize=8)
ax[1].set_xlabel('channels per cell', fontsize=8)
ax[0].set_ylabel('survival probability', fontsize=8)
ax[1].set_ylabel('survival probability', fontsize=8)
ax[0].set_title(r'slow shock (< 1.0 Hz)', backgroundcolor=color['pale_yellow'],
                y=1.01, fontsize=8)
ax[1].set_title(r'fast shock (≥ 1.0 Hz)', backgroundcolor=color['pale_yellow'],
                y=1.01, fontsize=8)


# Plot the data points.
_ = ax[0].plot(slow_exp['chan_per_cell'], slow_plot, 'k+', ms=1.5, alpha=0.5,
               label='individual cells')
_ = ax[1].plot(fast_exp['chan_per_cell'], fast_plot, 'k+', ms=1.5, alpha=0.5)


# Plot the regression curves.
_ = ax[0].plot(chan_range, slow_prob, color=color['red'],
               label='logistic regression')

_ = ax[1].plot(chan_range, fast_prob, color=color['blue'])


# Fill in the credible regions.
_ = ax[0].fill_between(chan_range, slow_cred_region[0, :], slow_cred_region[1, :],
                       color=color['light_red'], alpha=0.5)

_ = ax[0].plot(chan_range, slow_cred_region[0, :],
               color=color['red'], lw=0.75, alpha=0.8)
_ = ax[0].plot(chan_range, slow_cred_region[1, :],
               color=color['red'], lw=0.75, alpha=0.8)

_ = ax[1].fill_between(chan_range, fast_cred_region[0, :], fast_cred_region[1, :],
                       color=color['light_blue'], alpha=0.5)
_ = ax[1].plot(chan_range, fast_cred_region[0, :],
               color=color['blue'], lw=0.75, alpha=0.6)
_ = ax[1].plot(chan_range, fast_cred_region[1, :],
               color=color['blue'], lw=0.75, alpha=0.6)

for a in ax:
    a.tick_params(labelsize=8)
    a.set_xlim([0, 800])

ax[0].text(-0.17, 1.04, '(A)', fontsize=10, transform=ax[0].transAxes)
ax[1].text(-0.17, 1.04, '(B)', fontsize=10, transform=ax[1].transAxes)
plt.tight_layout()
plt.savefig('fig4.pdf', bbox_inches='tight')
plt.savefig('fig4.png', bbox_inches='tight')

# %%
shock_exp.loc[:, 'bin_number'] = 0
bins = np.arange(0, 800, 100)
dfs = []
for i in range(1, len(bins)):
    sub_set = shock_exp[(shock_exp['chan_per_cell'] >= bins[i - 1])
                        & (shock_exp['chan_per_cell'] < bins[i])]
    sub_set.loc[:, 'bin_number'] = i
    dfs.append(sub_set)

binned_df = pd.concat(dfs, ignore_index=True)

# %%
fig, ax = plt.subplots()
grouped = binned_df.groupby('bin_number')
for g, d in grouped:
    mean_chan = d['chan_per_cell'].mean()
    _grouped = d.groupby('flow_rate')
    surv_prob = []
    flow_rate = []
    for _g, _d in _grouped:
        _surv = _d['survival'].sum() / len(_d)
        surv_prob.append(_surv)
        flow_rate.append(_g)
    ax.plot(flow_rate, surv_prob, 's-', label=mean_chan)
# ax.set_xscale('log')
ax.legend()


# %%
fig, ax = plt.subplots(4, 3)
ax = ax.ravel()
for a in ax:
    a.set_ylim([-0.01, 1.01])
    a.set_xlim([0, 500])
    a.tick_params(labelsize=8)
sorted_df = binned_df.sort_values('flow_rate')
axes = {a: b for b, a in enumerate(sorted_df['flow_rate'].unique())}
grouped = sorted_df.groupby('flow_rate')

for g, d in grouped:
    surv = d['survival'].values.astype(
        int) - np.random.normal(loc=0, scale=0.01, size=len(d))
    ax[axes[g]].plot(d['chan_per_cell'], surv, 'k.')
    ax[axes[g]].set_title('{0} hz'.format(g), fontsize=8)

shock_exp[shock_exp['flow_rate'] == 0.02]
plt.tight_layout()
