# -*- coding: utf-8 -*-
# %%
import os
os.chdir('code/figs/')
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pystan
import datetime
sys.path.insert(0, '../../')
import mscl.plotting
import mscl.stats
import mscl.mcmc
colors = mscl.plotting.set_plotting_style()

FIG_NO = 'X_alternative_logits'

# Load the data
data = pd.read_csv('../../data/csv/mscl_survival_data.csv')

# Set the generic shock rate idx.
data['idx'] = 1
data.loc[data['shock_class'] == 'fast', 'idx'] = 2

data
#%%# Load the generic logistic regression stan model.
model = pystan.StanModel('../stan/generic_logistic_regression.stan')

# %% Use area as a predictor.
area_predictor_dict = {'J': 2, 'N': len(
    data), 'trial': data['idx'], 'output': data['survival'], 'predictor': data['area']}
area_samples = model.sampling(data=area_predictor_dict, iter=5000, chains=4)

# Compute the statistics.
area_df = mscl.mcmc.chains_to_dataframe(area_samples)
area_stats = mscl.mcmc.compute_statistics(area_df)

# %% Use shock rate as a predictor.
rate_predictor_dict = {'J': 1, 'N': len(data), 'trial': np.ones(
    len(data)).astype(int), 'predictor': data['flow_rate'], 'output': data['survival']}
rate_samples = model.sampling(data=rate_predictor_dict, iter=5000,  chains=4)
rate_samples
rate_df = mscl.mcmc.chains_to_dataframe(rate_samples)
rate_stats = mscl.mcmc.compute_statistics(rate_df)

# %% Use date as a predictor
# Reset the date to days.
date_list = np.sort(data['date'].unique())
date0 = datetime.datetime.strptime(str(date_list[0]), '%Y%m%d')
date_key = {date_list[0]: 0}
for i in range(1, len(date_list)):
    parsed_date = datetime.datetime.strptime(str(date_list[i]), '%Y%m%d')
    delta = parsed_date - date0
    date_key[date_list[i]] = delta.days

# Reassign the dates.
days = [date_key[i] for i in data['date']]
data.loc[:, 'delta_days'] = days

# Load the integer predictor model
date_model = pystan.StanModel(
    '../stan/integer_predictor_logistic_regression.stan')

#%% Sample
date_predictor_dict = {'J': 1, 'N': len(data), 'trial': np.zeros(len(data)).astype(
    int), 'output': data['survival'], 'predictor': data['delta_days']}
date_samples = date_model.sampling(data=date_predictor_dict, iter=100, chains=4)
date_samples

# %%
area_range = np.linspace(0, 50, 500)
shock_range = np.linspace(0, 2.5, 500)
predictor_range = [area_range, shock_range]
color_key = ['purple', 'blue']

# Generate the plot.
fig, ax = plt.subplots(2, 2, figsize=(6, 6), sharey=True)
ax[0, 0].set_title('slow shock (< 1.0 Hz)', fontsize=10,
                   backgroundcolor=colors['pale_yellow'], y=1.05)
ax[0, 1].set_title('fast shock ( $\geq$ 1.0 Hz)', fontsize=10,
                   backgroundcolor=colors['pale_yellow'], y=1.05)
for a in ax.ravel():
    a.set_ylim([0, 1])
for i in range(2):
    ax[0, i].set_xlabel('segmented  cell area [µm$^2$]')
    ax[1, i].set_xlabel('shock rate [Hz]')

samples = [area_stats, rate_stats]
for a in ax.ravel()[:2]:
    a.set_xlim([0, 50])

for a in ax.ravel()[2:]:
    a.set_xlim([0, 2.5])


for i in range(2):
    beta0_median = area_stats[area_stats['parameter']
                              == 'beta_0__{}'.format(i)]['median'].values[0]
    beta1_median = area_stats[area_stats['parameter']
                              == 'beta_1__{}'.format(i)]['median'].values[0]
    prediction = (1 + np.exp(-beta0_median - beta1_median * area_range))**-1

    # plot the cred regions.
    cred_region = np.zeros((2, len(area_range)))
    beta0_samples = area_df['beta_0__{}'.format(i)]
    beta1_samples = area_df['beta_1__{}'.format(i)]
    for j, a in enumerate(area_range):
        prob = (1 + np.exp(-beta0_samples - beta1_samples * a))**-1
        cred_region[:, j] = mscl.mcmc.compute_hpd(prob, 0.95)
    ax[0, i].plot(area_range, prediction, lw=2, color=colors[color_key[i]])
    ax[0, i].fill_between(area_range, cred_region[0, :], cred_region[1, :],
                          alpha=0.5, color=colors['light_' + color_key[i]])
    ax[0, i].plot(area_range, cred_region[1, :],
                  lw=1, color=colors[color_key[i]])
    ax[0, i].plot(area_range, cred_region[0, :],
                  lw=1, color=colors[color_key[i]])


# Plot the shock rate predictions.
beta0_median = rate_stats[rate_stats['parameter']
                          == 'beta_0']['median'].values[0]
beta1_median = rate_stats[rate_stats['parameter']
                          == 'beta_1']['median'].values[0]
beta0_samples = rate_df['beta_0']
beta1_samples = rate_df['beta_0']

prediction = (1 + np.exp(-beta0_median - beta1_median * shock_range))**-1
cred_region = np.zeros((2, len(shock_range)))
for i, r in enumerate(shock_range):
    prob = (1 + np.exp(-beta0_samples - beta1_samples * r))**-1
    cred_region[:, i] = mscl.mcmc.compute_hpd(prob, 0.95)

ax[1, 0].plot(shock_range, prediction, color='slategray', lw=2)
ax[1, 0].fill_between(shock_range, cred_region[0, :],
                      cred_region[1, :], color='slategray', alpha=0.3)
ax[1, 0].plot(shock_range, cred_region[0, :], color='slategray', lw=1)
ax[1, 0].plot(shock_range, cred_region[1, :], color='slategray', lw=1)
plt.tight_layout()

print('hello')
