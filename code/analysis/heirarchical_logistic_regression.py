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
STAN_MODEL = '../stan/heirarchical_logistic_regression.stan'
sm = stan.StanModel(STAN_MODEL)

# Load and isolate shock data
data = pd.read_csv('../../data/csv/mscl_survival_data.csv')
shock_data = data[data['experiment'] == 'shock'].copy()


#%% Group by flow class and sample.


# --------------------------------------------------------------------------
# WARNING: ADD TEMPORARY DATA FOR ZERO CHANNELS -- MUST BE FIXED:
# _________________________________________________________________________
trace_df = []
grouped = shock_data.groupby('shock_class')
for g, d in grouped:

    chan = list(d['effective_channels'].values)
    chan_err = list(d['effective_channel_err'].values)
    surv = list(d['survival'].values)

    for i in range(100):
        chan.append(0)
        chan_err.append(1E-9)
        if np.random.rand() < 1E-3:
            surv.append(1)
        else:
            surv.append(0)

    # Sample the model.
    data_dict = {'N': len(surv), 'predictor': chan,
                 'predictor_err': chan_err,
                 'output': surv}
    trace = sm.sampling(data=data_dict, iter=2000, chains=4)
    # Convert traces to dataframe.
    _df = mscl.mcmc.chains_to_dataframe(trace, var_names=['beta_0', 'beta_1'])
    _df.insert(0, 'shock_speed', g)
    trace_df.append(_df)

# Concatenate the two data frames.
traces = pd.concat(trace_df, ignore_index=True)


# Save both to disk.
traces.to_csv('../../data/csv/heirarchical_logistic_regression_traces.csv',
              index=False)
