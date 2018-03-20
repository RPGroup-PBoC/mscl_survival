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
data = pd.read_csv('../data/csv/compiled_data.csv')
shock_data = data[data['class'] == 'shock'].copy()
shock_data.loc[shock_data['flow_rate'] < 1.0, 'shock_speed'] = 'slow'
shock_data.loc[shock_data['flow_rate'] >= 1.0, 'shock_speed'] = 'fast'
sm = pystan.StanModel('../code/generic_logistic_regression.stan')

#%% Area as a predictor variable
dfs = []
grouped = shock_data.groupby('shock_speed')
for g, d in grouped:
    data_dict = {'N': len(d), 'predictor': d['area'],
                 'output': d['survival'].astype(int)}
    trace = sm.sampling(data=data_dict, iter=5000, chains=4)
    _df = mscl.mcmc.chains_to_dataframe(trace)
    _df.insert(0, 'shock_speed', g)
    dfs.append(_df)
traces = pd.concat(dfs, ignore_index=False)

traces.to_csv('../data/csv/logistic_regression_area_predictor.csv',
              index=False)

# %%
grouped = shock_data.groupby('shock_speed')
fig = plt.figure(figsize=(6.5, 3.5))
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
fig.text(0.5, 0.9, '(B)', fontsize=8)
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
plt.savefig('figS{}.png'.format(FIG_NO), bbox_inches='tight')
plt.savefig('figS{}.pdf'.format(FIG_NO), bbox_inches='tight')
