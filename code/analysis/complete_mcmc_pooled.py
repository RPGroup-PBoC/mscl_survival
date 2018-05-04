# -*- coding: utf-8 -*-
# %%
import sys
import numpy as np
import pandas as pd
import pystan
sys.path.insert(0, '../../')
import mscl.stats
import mscl.mcmc
import glob

# Define the parameters from Bialecka-Fornal et al.  2012
CHANNEL_NUM = 340
CHANNEL_SEM = 64

# load the data sets.
files = glob.glob('../processing/*mlg910*/output/*.csv')
mlg910 = pd.concat([pd.read_csv(f, comment='#')
                    for f in files], ignore_index=True)
# Adjust the area for the 2018 measurement.
mlg910.loc[mlg910['date'] == 20180410,
           'area'] = mlg910[mlg910['date'] == 20180410]['area'] * 0.5

# Adjust the replicate number for the July 2017 run.
mlg910.loc[mlg910['date'] == 20170721,
           'replicate_number'] = mlg910['replicate_number'].max() + 1
mlg910.loc[:, 'survival'] = -1
shock_files = glob.glob('../processing/*shock*/output/*.csv')
shock_data = pd.concat([pd.read_csv(f, comment='#') for f in shock_files])
shock_data['replicate_number'] = 0
data = pd.concat([mlg910, shock_data], ignore_index=True)
data.dropna(axis=1, inplace=True)

# Insert the shock rate identifier.
data.loc[:, 'idx'] = 0

# Rescale the intensities to meaningful values.
max_exp = data['exposure_ms'].max()
data['scaled_intensity'] = (
    data['intensity'] - data['mean_bg']) * max_exp / data['exposure_ms']

data = data[data['scaled_intensity'] >= 0].copy()

# Load the stan model.
model = pystan.StanModel('../stan/complete_analysis.stan')


# Assemble the data dictionary.
data_dict = dict(J1=6, J2=1, N1=len(data[data['rbs'] == 'mlg910']), N2=len(data[data['rbs'] != 'mlg910']),
                 repl=data[data['replicate_number'] > 0]['replicate_number'].values.astype(int),
                 shock=np.ones(len(data[data['rbs'] != 'mlg910'])).astype(int), 
                 A=data[data['rbs']=='mlg910']['area'],
                 I_A_sc=data[data['rbs'] == 'mlg910']['scaled_intensity'],
                 I_A_sr=data[data['rbs'] != 'mlg910']['scaled_intensity'],
                 ntot=CHANNEL_NUM, sigma_ntot=CHANNEL_SEM,
                 survival=data[data['rbs'] != 'mlg910']['survival'].values.astype(int))
#%% Sample
print('sampling...')
samples = model.sampling(data=data_dict, iter=8000, chains=4, n_jobs=-1)
print('finished!')

# %% Convert to a DataFrame and and save.
sample_df = mscl.mcmc.chains_to_dataframe(samples)
stats = mscl.mcmc.compute_statistics(sample_df)
sample_df.to_csv('../../data/csv/pooled_complete_mcmc_traces.csv', index=False)
stats.to_csv('../../data/csv/pooled_complete_mcmc_stats.csv', index=False)

# Piece together the final data frame
shock_data = data[data['rbs'] != 'mlg910'].copy()
cal_data = data[data['rbs']=='mlg910'].copy()
shock_data.drop(axis=1, labels=['exposure_ms', 'mean_bg', 'idx', 'intensity', 'replicate_number'], inplace=True)

for j in range(2):
	beta_stats = stats[stats['parameter'] == 'beta_{}'.format(j)].values[0][1:]
	shock_data.loc[:, 'beta_{}_median'.format(j)] = beta_stats[0] 
	shock_data.loc[:, 'beta_{}_min'.format(j)] = beta_stats[1] 
	shock_data.loc[:, 'beta_{}_max'.format(j)] = beta_stats[2] 

# Include the channel copy number info.
median_n = [stats[stats['parameter']=='n_sr__{}'.format(i)]['median'].values[0] for i in range(len(shock_data))]
min_n = [stats[stats['parameter']=='n_sr__{}'.format(i)]['hpd_min'].values[0] for i in range(len(shock_data))]
max_n = [stats[stats['parameter']=='n_sr__{}'.format(i)]['hpd_max'].values[0] for i in range(len(shock_data))]
shock_data.loc[:, 'effective_channels'] = median_n
shock_data.loc[:, 'minimum_channels'] = min_n
shock_data.loc[:, 'maximum_channels'] = max_n

shock_data.loc[:, 'calibration_factor'] = stats[stats['parameter']=='hyper_alpha_mu']['median'].values[0]
shock_data.loc[:, 'minimum_calibration_factor'] = stats[stats['parameter']=='hyper_alpha_mu']['hpd_min'].values[0]
shock_data.loc[:, 'maximum_calibration_factor'] = stats[stats['parameter']=='hyper_alpha_mu']['hpd_max'].values[0]

# Save the final dataframe.
shock_data.to_csv('../../data/csv/pooled_mscl_survival_data.csv', index=False)
