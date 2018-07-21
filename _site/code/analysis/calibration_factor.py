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
model = pystan.StanModel('../stan/hierarchical_calibration_factor.stan')

# Set the data dictionary.
data_dict = {'J': int(mlg910['replicate_number'].max()), 'N': len(mlg910), 'rep': mlg910['replicate_number'].values.astype(int),
             'A': mlg910['area'], 'I_A': mlg910['scaled_intensity'], 'ntot': CHANNEL_NUM, 'sigma_ntot': CHANNEL_SEM}

samples = model.sampling(data=data_dict, iter=100000, chains=4, thin=10)
df = mscl.mcmc.chains_to_dataframe(samples)

# Compute the median and HPDs.
alpha_mu = np.median(df['hyper_alpha_mu'])
alpha_hpd = mscl.mcmc.compute_hpd(df['hyper_alpha_mu'], 0.95)
A_mu = np.median(df['hyper_A_mu'])
A_hpd = mscl.mcmc.compute_hpd(df['hyper_A_mu'], 0.95)


# Plot the distribution of area and alpha.
fig, ax = plt.subplots(1, 1, figsize=(5, 5))
_ = sns.kdeplot(df['hyper_alpha_mu'], df['hyper_A_mu'],
                shade=True, cmap='Blues', ax=ax)
ax.set_xlabel(r'$\tilde{\alpha}$ [a.u. / channel]')
ax.set_ylabel(r'$\tilde{A}$ [$\mu$m$^2$]')
ax.set_ylim([1, 8])
ax.set_xlim([750, 2250])
_ = ax.plot(alpha_mu, A_mu, 'o', color=colors['red'])
_ = ax.hlines(A_mu, alpha_hpd[0], alpha_hpd[1], color=colors['red'])
_ = ax.vlines(alpha_mu, A_hpd[0], A_hpd[1], color=colors['red'])
ax.set_title(r'$\tilde{\alpha} = %d^{+%d}_{-%d}$  $\tilde{A} = %.2f^{+%.2f}_{-%.2f}$' % (alpha_mu,
            alpha_hpd[1] - alpha_mu, alpha_mu - alpha_hpd[0], A_mu, A_hpd[1] - A_mu, A_mu - A_hpd[0]), 
            backgroundcolor=colors['pale_yellow'])
plt.savefig('../../figs/hyperparameter_posterior.png', bbox_inches='tight')
