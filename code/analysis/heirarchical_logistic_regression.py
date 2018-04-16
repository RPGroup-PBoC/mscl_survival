# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pystan as stan
import sys
sys.path.insert(0, '../../')
import mscl.plotting
import mscl.mcmc
import imp
# imp.reload(mscl.mcmc)
colors = mscl.plotting.set_plotting_style()

<<<<<<< HEAD
ALPHA = 1434
ALPHA_SIG = 504
A = 4.85
A_SIG = 1.27 
=======
>>>>>>> 709a6946ca985c61e256071636766d53e0e2ee5d
# Define and load stan model.
STAN_MODEL = '../stan/hierarchical_logistic_regression.stan'
sm = stan.StanModel(STAN_MODEL)

<<<<<<< HEAD
# Load the data and assemble the dictionary.  
data = pd.read_csv('../../data/csv/compiled_data.csv')

# Add the identifier column.
data.loc[:, 'idx'] = 1
data.loc[data['shock_class']=='fast', 'idx'] = 2

data_dict = {'J': 2, 'N':len(data), 'trial':data['idx'], 'I_A':data['rescaled_intensity'],
            'alpha': ALPHA, 'alpha_sigma': ALPHA_SIG, 'A':A, 'A_sigma':A_SIG, 
            'survival': data['survival'].values.astype(int)}
samples = sm.sampling(data=data_dict, iter=100000, chains=4, thin=100)

df = mscl.mcmc.chains_to_dataframe(samples)
stats = mscl.mcmc.compute_statistics(df)

# Assemble the complete dataframe.
median_channels = [stats[stats['parameter']=='n__{}'.format(i)]['median'].values[0] for i in range(len(data))] 
min_channels = [stats[stats['parameter']=='n__{}'.format(i)]['hpd_min'].values[0] for i in range(len(data))] 
max_channels = [stats[stats['parameter']=='n__{}'.format(i)]['hpd_max'].values[0] for i in range(len(data))] 
data['effeective_channels'] = median_channels
data['min_channels'] = min_channels
data['max_channels'] = max_channels
data.to_csv('../../data/csv/mscl_survival_data.csv', index=False)

# Save the samples from the logistic regression.
logit_df = df[['beta_0__0', 'beta_0__1', 'beta_1__0', 'beta_1__1']]
logit_df.columns = ['beta_0__slow', 'beta_0__fast', 'beta_1__slow', 'beta_1__fast']
logit_df.to_csv('../../data/csv/logistic_regression_traces.csv', index=False)
=======
# Load and isolate shock data
shock_data = pd.read_csv('../../data/csv/compiled_data.csv')
# shock_data = data[data['experiment'] == 'shock'].copy()

#%% Perform regression
data_dict = {'N': len(shock_data), 'predictor': shock_data['effective_channels'],
             'output': shock_data['survival'].astype(int)}
trace = sm.sampling(data=data_dict, iter=2000, chains=4)
pooled_df = mscl.mcmc.chains_to_dataframe(
    trace, var_names=['beta_0', 'beta_1'])
pooled_df.to_csv('../../data/csv/pooled_generic_logistic_regression.csv')

#%% Perform regression
data_dict = {'N': len(shock_data), 'predictor': np.log(shock_data['effective_channels']),
             'output': shock_data['survival'].astype(int)}
trace = sm.sampling(data=data_dict, iter=2000, chains=4)
pooled_df = mscl.mcmc.chains_to_dataframe(
    trace, var_names=['beta_0', 'beta_1'])
pooled_df.to_csv('../../data/csv/pooled_generic_logistic_regression_log.csv')

#%% Group by flow class and sample.
# --------------------------------------------------------------------------
# WARNING: ADD TEMPORARY DATA FOR ZERO CHANNELS -- MUST BE FIXED:
# _________________________________________________________________________
trace_df = []
grouped = shock_data.groupby('shock_class')
for g, d in grouped:

    data_dict = {'N': len(d), 'predictor': np.log(d['effective_channels']),
                 'output': d['survival'].astype(int)}
    trace = sm.sampling(data=data_dict, iter=2000, chains=4)

    # Convert traces to dataframe.
    _df = mscl.mcmc.chains_to_dataframe(trace, var_names=['beta_0', 'beta_1'])
    _df.insert(0, 'shock_speed', g)
    trace_df.append(_df)

traces = pd.concat(trace_df, ignore_index=True)
# Save both to disk.
traces.to_csv('../../data/csv/heirarchical_logistic_regression_traces.csv',
              index=False)
>>>>>>> 709a6946ca985c61e256071636766d53e0e2ee5d
