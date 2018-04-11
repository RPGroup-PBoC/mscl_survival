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

# Define and load stan model.
STAN_MODEL = '../stan/generic_logistic_regression.stan'
sm = stan.StanModel(STAN_MODEL)

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
