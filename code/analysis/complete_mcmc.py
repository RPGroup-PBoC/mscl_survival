# -*- coding: utf-8 -*-
#%%
import numpy as np
import pandas as pd
import pystan
import sys
import glob
sys.path.insert(0, '../../')
import mscl.stats
import mscl.mcmc

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
mlg910.loc[:, 'shock_class'] = 'none'
mlg910.loc[:, 'survival'] = -1
shock_files = glob.glob('../processing/*shock*/output/*.csv')
shock_data = pd.concat([pd.read_csv(f, comment='#') for f in shock_files])
shock_data['replicate_number'] = 0
shock_data.loc[:, 'shock_class'] = 'slow'
shock_data.loc[shock_data['flow_rate'] >= 1.0, 'shock_class'] = 'fast'
data = pd.concat([mlg910, shock_data], ignore_index=True)
data.dropna(axis=1, inplace=True)
# Insert the shock rate identifier.
data.loc[:, 'idx'] = 0
data.loc[data['shock_class'] == 'slow', 'idx'] = 1
data.loc[data['shock_class'] == 'fast', 'idx'] = 2

# Rescale the intensities to meaningful values.
max_exp = data['exposure_ms'].max()
data['scaled_intensity'] = (
    data['intensity'] - data['mean_bg']) * max_exp / data['exposure_ms']

data = data[data['scaled_intensity'] >= 0].copy()

# Load the stan model.
model = pystan.StanModel('../stan/complete_analysis.stan')


# Assemble the data dictionary.
data_dict = dict(J1=6, J2=2, N1=len(data[data['rbs'] == 'mlg910']), N2=len(data[data['rbs'] != 'mlg910']),
                 repl=data[data['replicate_number'] >
                           0]['replicate_number'].values.astype(int),
                 shock=data[data['idx'] >
                            0]['idx'], A=data[data['rbs'] == 'mlg910']['area'],
                 I_A_sc=data[data['rbs'] == 'mlg910']['scaled_intensity'],
                 I_A_sr=data[data['rbs'] != 'mlg910']['scaled_intensity'],
                 ntot=CHANNEL_NUM, sigma_ntot=CHANNEL_SEM,
                 survival=data[data['rbs'] != 'mlg910']['survival'].values.astype(int))

#%% Sample
print('sampling...')
samples = model.sampling(data=data_dict, iter=10000, chains=3)
print('finished!')

#%% Compute the statistics of each parameter.
sample_df = mscl.mcmc.chains_to_dataframe(samples)
stats = mscl.mcmc.compute_statistics(sample_df)
sample_df.to_csv('../../data/csv/complete_mcmc_traces.csv', index=False)
stats.to_csv('../../data/csv/complete_mcmc_stats.csv', index=False)

# Piece together the final data frame
shock_data = data[data['rbs'] != 'mlg910'].copy()
cal_data = data[data['rbs']=='mlg910'].copy()
shock_data.drop(axis=1, labels=['exposure_ms', 'mean_bg', 'idx', 'intensity', 'replicate_number'], inplace=True)
cal_data.drop(axis=1, labels=['exposure_ms', 'mean_bg', 'idx', 'intensity', 'replicate_number', 'shock_class'], inplace=True)

rate_keys = {'slow':0, 'fast':1}
for i, r in enumerate(rate_keys.keys()):
    for j in range(2):
        beta_stats = stats[stats['parameter'] == 'beta_{}__{}'.format(j, rate_keys[r])].values[0][1:]
        shock_data.loc[shock_data['shock_class']==r, 'beta_{}_median'.format(j)] = beta_stats[0] 
        shock_data.loc[shock_data['shock_class']==r, 'beta_{}_min'.format(j)] = beta_stats[1] 
        shock_data.loc[shock_data['shock_class']==r, 'beta_{}_max'.format(j)] = beta_stats[2] 

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

# median_n = [stats[stats['parameter']=='n_sc__{}'.format(i)]['median'].values[0] for i in range(6)]
# min_n = [stats[stats['parameter']=='n_sc__{}'.format(i)]['hpd_min'].values[0] for i in range(6)]
# max_n = [stats[stats['parameter']=='n_sc__{}'.format(i)]['hpd_max'].values[0] for i in range(6)]
# cal_data.loc[:, 'effective_channels'] = median_n
# cal_data.loc[:, 'minimum_channels'] = min_n
# cal_data.loc[:, 'maximum_channels'] = max_n

cal_data.loc[:, 'calibration_factor'] = stats[stats['parameter']=='hyper_alpha_mu']['median'].values[0]
cal_data.loc[:, 'minimum_calibration_factor'] = stats[stats['parameter']=='hyper_alpha_mu']['hpd_min'].values[0]
cal_data.loc[:, 'maximum_calibration_factor'] = stats[stats['parameter']=='hyper_alpha_mu']['hpd_max'].values[0]
cal_data.drop(axis=1, labels=['exposure_ms', 'mean_bg', 'idx', 'intensity', 'replicate_number', 'shock_class'], inplace=True)
#
# Save the final dataframe.
shock_data.to_csv('../../data/csv/mscl_survival_data.csv', index=False)
cal_data.to_csv('../../data/csv/standard_candle_data.csv', index=False)
