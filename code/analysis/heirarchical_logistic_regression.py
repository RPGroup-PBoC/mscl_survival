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

colors = mscl.plotting.set_plotting_style()
%matplotlib inline
# Define and load stan model.
STAN_MODEL = '../heirarchical_logistic_regression.stan'
sm = mscl.mcmc.load_StanModel(STAN_MODEL)

# Load and isolate shock data
data = pd.read_csv('../../data/csv/compiled_data.csv')
shock_data = data[data['class'] == 'shock'].copy()

# Separate by shock speed.
shock_data.loc[shock_data['flow_rate'] >= 1.0, 'shock_speed'] = 'fast'
shock_data.loc[shock_data['flow_rate'] < 1.0, 'shock_speed'] = 'slow'

#%% Group by flow class and sample.
trace_df = []
grouped = shock_data.groupby('shock_speed')
for g, d in grouped:
    # Sample the model.
    data_dict = {'N': len(d), 'num_channel': d['chan_per_cell'],
                 'channel_err': 0.15 * d['chan_per_cell'],
                 'survival': d['survival'].astype(int)}
    trace = sm.sampling(data=data_dict, iter=2000, chains=4)
    # Convert traces to dataframe.
    _df = mscl.mcmc.chains_to_dataframe(trace, var_names=['alpha', 'beta'])
    _df.insert(0, 'shock_speed', g)
    trace_df.append(_df)

# Concatenate the two data frames.
traces = pd.concat(trace_df, ignore_index=True)


# Save both to disk.
traces.to_csv('../../data/csv/heirarchical_logistic_regression_traces.csv',
              index=False)
