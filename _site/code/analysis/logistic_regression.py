# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import pystan
sys.path.insert(0, '../../')
import mscl.mcmc
import mscl.stats

# Load the data
data = pd.read_csv('../../data/csv/compiled_data.csv')
data.loc[:, 'idx'] = 1
data.loc[data['shock_class'] == 'fast', 'idx'] = 2
data = data[data['effective_channels'] > 0]

# Compile the stan model.
model = pystan.StanModel(
    '../stan/heirarchical_logistic_regression_error_propagation.stan')

# Assemble the data dictionary.
data_dict = {'J': 2, 'N': len(data), 'trial': data['idx'], 'n_c': data['effective_channels'],
             'sigma_c': data['effective_channels_err'], 'survival': data['survival'].values.astype(int)}

# Sample
samples = model.sampling(data=data_dict, iter=5000, chains=4)
sample_df = mscl.mcmc.chains_to_dataframe(samples)
idx_dict = {1:'slow', 2:'fast'}
new_cols = {'beta_{}__{}'.format(i, j):'beta_{}_{}'.format(i, idx_dict[j+1]) for i in range(2) for j in range(2)}
sample_df.rename(columns=new_cols, inplace=True)
sample_df.to_csv('../../data/csv/logistic_regression_traces.csv', index=False)