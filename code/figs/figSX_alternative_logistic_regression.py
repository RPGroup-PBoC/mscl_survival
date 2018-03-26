# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import pystan
import pandas as pd
import sys
sys.path.insert(0, '../')
import mscl.stats
import mscl.mcmc
import mscl.plotting
colors = mscl.plotting.set_plotting_style()

FIG_NO = 'X3'

# --- Area Predictor -------------
data = pd.read_csv('../../data/csv/mscl_survival_data.csv')
shock_data = data[data['experiment'] == 'shock'].copy()
sm = pystan.StanModel('../stan/generic_logistic_regression.stan')


#%% Area as a predictor variable
dfs = []
grouped = shock_data.groupby('shock_class')
for g, d in grouped:
    data_dict = {'N': len(d), 'predictor': d['area'],
                 'output': d['survival'].astype(int)}
    trace = sm.sampling(data=data_dict, iter=5000, chains=4)
    _df = mscl.mcmc.chains_to_dataframe(trace)
    _df.insert(0, 'shock_speed', g)
    dfs.append(_df)
traces = pd.concat(dfs, ignore_index=False)

traces.to_csv('../../data/csv/logistic_regression_area_predictor.csv',
              index=False)

# %%
grouped = shock_data.groupby('shock_class')
fig = plt.figure(figsize=(4, 2.2))
gs = gridspec.GridSpec(3, 2, height_ratios=[0.5, 8, 0.5])
ax = [plt.subplot(gs[i, j]) for i in range(3) for j in range(2)]
axes = {'slow': [ax[0], ax[4], ax[2]], 'fast': [ax[1], ax[5], ax[3]]}
color_dict = {'slow': [colors['red'], colors['light_red']],
              'fast': [colors['blue'], colors['light_blue']]}


ax[0].set_title('slow shock (< 1.0 Hz)', backgroundcolor=colors['pale_yellow'],
                fontsize=8)
ax[1].set_title('fast shock (≥ 1.0 Hz)', backgroundcolor=colors['pale_yellow'],
                fontsize=8)
fig.text(0.05, 0.9, '(A)', fontsize=8)
# Set the special axes requirements for the rug plots
for i in (0, 1):
    ax[i].set_axis_off()

for i in (4, 5):
    ax[i].set_frame_on(False)
    ax[i].set_yticks([])
    ax[i].set_xlabel('area [µm$^2$]', fontsize=8)

for a in ax:
    a.tick_params(labelsize=8)

ax[2].set_ylabel('survival probability', fontsize=8)
ax[3].set_ylabel('survival probability', fontsize=8)
ax[2].set_xticklabels([])
ax[3].set_xticklabels([])

# Plot the survival and death cells on the appropriate rug plots
area_range = np.linspace(0, 40, 1000)
for g, d in grouped:
    # Generate the rug plots
    survival = d[d['survival'] == True]
    death = d[d['survival'] == False]
    _y_surv = np.random.normal(loc=0, scale=0.1, size=len(survival))
    _y_death = np.random.normal(loc=0, scale=0.1, size=len(death))
    surv_ax, death_ax, plot_ax = axes[g]
    surv_ax.plot(survival['area'], _y_surv, 'k.', ms=1.5, alpha=0.2)
    death_ax.plot(death['area'], _y_death, 'k.', ms=1.5, alpha=0.2)

grouped = traces.groupby('shock_speed')
for g, d in grouped:
    _, _, plot_ax = axes[g]
    # Compute the survival probability curve.
    beta_0, beta_1 = np.median(d['beta_0']), np.median(d['beta_1'])
    p_survival = (1 + np.exp(-beta_0 - beta_1 * area_range))**-1

    # Compute the credible region.
    cred_region = np.zeros((2, len(area_range)))
    for i, a in enumerate(area_range):
        _p_surv = (1 + np.exp(-d['beta_0'] - a * d['beta_1']))**-1
        cred_region[:, i] = mscl.mcmc.compute_hpd(_p_surv, 0.95)

    _ = plot_ax.plot(area_range, p_survival, '-', color=color_dict[g][0])
    _ = plot_ax.fill_between(area_range, cred_region[0, :], cred_region[1, :],
                             color=color_dict[g][1], alpha=0.5)
    _ = plot_ax.plot(area_range, cred_region[0, :], color=color_dict[g][0],
                     lw=1)
    _ = plot_ax.plot(area_range, cred_region[1, :], color=color_dict[g][0],
                     lw=1)

for a in ax:
    a.set_xlim([0, 40])
plt.subplots_adjust(hspace=0.01, wspace=0.3)
plt.savefig('../../figs/figS{}_area.pdf'.format(FIG_NO), bbox_inches='tight')
# plt.savefig('../figs/figS{}.pdf'.format(FIG_NO), bbox_inches='tight')
# %%
CAL = 4258
#%% Area as a predictor variable
dfs = []
grouped = shock_data.groupby('shock_class')
for g, d in grouped:
    data_dict = {'N': len(d), 'predictor': d['area'] * d['scaled_intensity'] / CAL,
                 'output': d['survival'].astype(int)}
    trace = sm.sampling(data=data_dict, iter=5000, chains=4)
    _df = mscl.mcmc.chains_to_dataframe(trace)
    _df.insert(0, 'shock_speed', g)
    dfs.append(_df)
traces = pd.concat(dfs, ignore_index=False)

# traces.to_csv('../../data/csv/logistic_regression_area_predictor.csv',
# index=False)

# %%
grouped = shock_data.groupby('shock_class')
fig = plt.figure(figsize=(4, 2.2))
gs = gridspec.GridSpec(3, 2, height_ratios=[0.5, 8, 0.5])
ax = [plt.subplot(gs[i, j]) for i in range(3) for j in range(2)]
axes = {'slow': [ax[0], ax[4], ax[2]], 'fast': [ax[1], ax[5], ax[3]]}
color_dict = {'slow': [colors['red'], colors['light_red']],
              'fast': [colors['blue'], colors['light_blue']]}


ax[0].set_title('slow shock (< 1.0 Hz)', backgroundcolor=colors['pale_yellow'],
                fontsize=8)
ax[1].set_title('fast shock (≥ 1.0 Hz)', backgroundcolor=colors['pale_yellow'],
                fontsize=8)
fig.text(0.05, 0.9, '(B)', fontsize=8)

# Set the special axes requirements for the rug plots
for i in (0, 1):
    ax[i].set_axis_off()

for i in (4, 5):
    ax[i].set_frame_on(False)
    ax[i].set_yticks([])
    ax[i].set_xlabel('integrated channel\ncopy number', fontsize=8)

for a in ax:
    a.tick_params(labelsize=8)

ax[2].set_ylabel('survival probability', fontsize=8)
ax[3].set_ylabel('survival probability', fontsize=8)
ax[2].set_xticklabels([])
ax[3].set_xticklabels([])

# Plot the survival and death cells on the appropriate rug plots
chan_range = np.linspace(0, 1500, 1000)
for g, d in grouped:
    # Generate the rug plots
    survival = d[d['survival'] == True]
    death = d[d['survival'] == False]
    _y_surv = np.random.normal(loc=0, scale=0.1, size=len(survival))
    _y_death = np.random.normal(loc=0, scale=0.1, size=len(death))
    surv_ax, death_ax, plot_ax = axes[g]
    surv_ax.plot(survival['area'] * survival['scaled_intensity'] /
                 CAL, _y_surv, 'k.', ms=1.5, alpha=0.2)
    death_ax.plot(death['area'] * death['scaled_intensity'] /
                  CAL, _y_death, 'k.', ms=1.5, alpha=0.2)

grouped = traces.groupby('shock_speed')
for g, d in grouped:
    _, _, plot_ax = axes[g]
    # Compute the survival probability curve.
    beta_0, beta_1 = np.median(d['beta_0']), np.median(d['beta_1'])
    p_survival = (1 + np.exp(-beta_0 - beta_1 * area_range))**-1

    # Compute the credible region.
    cred_region = np.zeros((2, len(area_range)))
    for i, a in enumerate(chan_range):
        _p_surv = (1 + np.exp(-d['beta_0'] - a * d['beta_1']))**-1
        cred_region[:, i] = mscl.mcmc.compute_hpd(_p_surv, 0.95)

    _ = plot_ax.plot(area_range, p_survival, '-', color=color_dict[g][0])
    _ = plot_ax.fill_between(chan_range, cred_region[0, :], cred_region[1, :],
                             color=color_dict[g][1], alpha=0.5)
    _ = plot_ax.plot(chan_range, cred_region[0, :], color=color_dict[g][0],
                     lw=1)
    _ = plot_ax.plot(chan_range, cred_region[1, :], color=color_dict[g][0],
                     lw=1)
for a in ax:
    a.set_xlim([0, 1500])
plt.subplots_adjust(hspace=0.01, wspace=0.3)
plt.savefig('../../figs/figS{}_total_channels.pdf'.format(FIG_NO),
            bbox_inches='tight')
# plt.savefig('../fi

# %%
sm = pystan.StanModel('../stan/bivariate_logistic_regression.stan')

#%% -- Bivariate regression --------------
# Compile the stan model.
sm = pystan.StanModel('../stan/bivariate_logistic_regression.stan')
data_dict = {'N': len(shock_data), 'predictor_1': shock_data['effective_channels'],
             'predictor_2': shock_data['flow_rate'],
             'output': shock_data['survival'].astype(int)}
# Sample and convert the traces to a dataframe
bivariate_traces = sm.sampling(data=data_dict, iter=5000, chains=4)
bivariate_df = mscl.mcmc.chains_to_dataframe(bivariate_traces)

# Compute the means for all coefficients.
_, beta_2, beta_1, beta_0 = np.mean(bivariate_df).values


#%% Set up the grid
chan_range = np.linspace(0, 800, 500)
shock_range = np.linspace(0, 3, 500)
X, Y = np.meshgrid(chan_range, shock_range)
p_survival = (1 + np.exp(-(beta_0 + beta_1 * X + beta_2 * Y)))**-1

# plt.imshow(p_survival)
fig, ax = plt.subplots(1, 1, figsize=(2.2, 2.2))
contour = ax.contourf(X, Y, p_survival, cmap='bone')
ax.set_yscale('log')

survs = shock_data[shock_data['survival'] == True]
deaths = shock_data[shock_data['survival'] == False]
# plt.clabel(contour, inline=1, fontsize=8)
ax.plot(survs['effective_channels'], survs['flow_rate'], '.', color=colors['green'],
        alpha=0.5, label='survival', ms=2)
ax.plot(deaths['effective_channels'], deaths['flow_rate'], '.', color=colors['light_purple'],
        alpha=0.5, label='death', ms=2)
cbar = plt.colorbar(mappable=contour, ax=ax)
ax.set_xlabel('effecive channel number', fontsize=8)
ax.set_ylabel('flow rate [Hz]', fontsize=8)
ax.xaxis.set_tick_params(labelsize=8)
ax.yaxis.set_tick_params(labelsize=8)
ax.legend(fontsize=8, bbox_to_anchor=(1.1, 1.2), ncol=2)
cbar.ax.yaxis.set_tick_params(labelsize=8)
cbar.ax.set_ylim([0, 1])
cbar.ax.set_ylabel('survival probability', fontsize=8)
plt.savefig('../../figs/figS{}_bivariate_logistic.pdf'.format(FIG_NO),
            bbox_inches='tight')
