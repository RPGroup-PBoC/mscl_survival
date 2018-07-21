# -*- coding: utf-8 -*-
# %%
import sys
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pystan
sys.path.insert(0, '../../')
import mscl.mcmc
import mscl.stats
import mscl.plotting
import seaborn as sns
colors = mscl.plotting.set_plotting_style()

# Define the parameters from Bialecka-Fornal et al.  2012
CHANNEL_NUM = 340
CHANNEL_SEM = 64

# load the data sets.
files = glob.glob('../processing/*mlg910*/output/*.csv')
mlg910 = pd.concat([pd.read_csv(f, comment='#')
                    for f in files], ignore_index=True)

# Adjust the area for the 2018 measurement.
mlg910.loc[:, 'area'] = mlg910['area'] / 2
# Ensure the exposure is set correctly.
max_exp = mlg910['exposure_ms'].max()
mlg910.loc[:, 'scaled_intensity'] = (
    mlg910['intensity'] - mlg910['mean_bg']) * max_exp / mlg910['exposure_ms']

mlg910 = mlg910[mlg910['scaled_intensity'] >= 0]
# Adjust the replicate number for the July 2017 run.
mlg910.loc[mlg910['date'] == 20170721,
           'replicate_number'] = mlg910['replicate_number'].max() + 1

# Load the stan model.
model = pystan.StanModel('../stan/heirarchical_standard_candle.stan')

# Assemble the data dictionary.
data_dict = {'J': int(mlg910['replicate_number'].max()), 'N': len(mlg910),
             'replicate': mlg910['replicate_number'].values.astype(int),
             'area_measurement': mlg910['area'], 'expression_measurement': mlg910['scaled_intensity']}
sm = model.sampling(data=data_dict, iter=5000, chains=4)

# save to a dataframe and compute the parameters
extract = sm.extract()
area_mu = extract['hyper_mu_area']
expression_mu = extract['hyper_mu_expression']
area_hpd = mscl.mcmc.compute_hpd(area_mu, 0.95)
expression_hpd = mscl.mcmc.compute_hpd(expression_mu, 0.95)

# Plot the corner plot.
fig, ax = plt.subplots(2, 2, figsize=(6, 6))
axes = [ax[0, 0], ax[1, 0], ax[1, 1]]
ax[0, 1].axis('off')

_ = axes[0].hist(area_mu, bins=100, normed=True, color=colors['red'])
_ = axes[2].hist(expression_mu, bins=100, normed=True, color=colors['red'])
_ = sns.kdeplot(area_mu, expression_mu, cbar=False,
                shade=True, cmap='inferno', ax=axes[1])

axes[0].set_yticks([])
axes[2].set_yticks([])

axes[0].set_title('area', backgroundcolor=colors['pale_yellow'], y=1.05)
axes[1].set_xlabel('area (µm$^2$')
axes[1].set_ylabel('expression [a.u.]')
axes[2].set_title('expression', backgroundcolor=colors['pale_yellow'], y=1.05)
plt.tight_layout()
plt.savefig('../../figs/standard_candle_mcmc.pdf', bbox_inches='tight')


# Compute parameters for the standard candle.
au_per_channel = np.median(area_mu) * np.median(expression_mu) / CHANNEL_NUM
au_per_channel_err = au_per_channel * ((CHANNEL_SEM / CHANNEL_NUM) +
                                       (np.std(expression_mu) / np.median(expression_mu)) + (np.std(area_mu) / np.median(area_mu)))

# Save the results in an easy to parse data frame
calibration_factor = pd.DataFrame(np.array(
    [au_per_channel, au_per_channel_err, np.median(area_mu), np.std(area_mu), max_exp])).T
calibration_factor.columns = ['au_per_channel',
                              'au_per_channel_err', 'mean_area', 'mean_area_err', 'exposure_ms']

calibration_factor.to_csv('../../data/csv/calibration_factor.csv', index=False)
