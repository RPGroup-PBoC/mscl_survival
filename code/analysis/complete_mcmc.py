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
mlg910.loc[mlg910['date'] == 20180410, 'area'] = mlg910[mlg910['date']==20180410]['area'] * 0.5

# Adjust the replicate number for the July 2017 run.
mlg910.loc[mlg910['date'] == 20170721,
           'replicate_number'] = mlg910['replicate_number'].max() + 1
mlg910.loc[:, 'shock_class'] = 'none'
mlg910.loc[:, 'survival'] = -1
shock_data = pd.read_csv('../../data/csv/compiled_data.csv')
shock_data['replicate_number'] = 0
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
data_dict = dict(J1=6, J2=2, N1=len(data[data['rbs']=='mlg910']), N2=len(data[data['rbs'] != 'mlg910']),
                 repl=data[data['replicate_number'] > 0]['replicate_number'].values.astype(int),
                 shock=data[data['idx'] > 0]['idx'], A=data[data['rbs'] == 'mlg910']['area'], 
                 I_A_sc=data[data['rbs'] == 'mlg910']['scaled_intensity'],
                 I_A_sr=data[data['rbs'] != 'mlg910']['scaled_intensity'],
                 ntot=CHANNEL_NUM, sigma_ntot=CHANNEL_SEM, 
                 survival=data[data['rbs']!= 'mlg910']['survival'].values.astype(int))


samples = model.sampling(data=data_dict, iter=5000, thin=10, chains=2)