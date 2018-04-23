# %%
# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.stats
import seaborn as sns
sys.path.insert(0, '../../')
import mscl.stats
import mscl.plotting
colors = mscl.plotting.set_plotting_style()


# Load the data sets.
shock_data = pd.read_csv('../../data/csv/mscl_survival_data.csv')
cal_data = pd.read_csv('../../data/csv/standard_candle_data.csv')
data = pd.concat([shock_data, cal_data])

# Compute the rescaled intensity relative to mlg910.
mean_mlg910 = np.mean(cal_data['scaled_intensity'])
data['relative_intensity'] = data['scaled_intensity'] / mean_mlg910

# Group by rbs.
grouped = data.groupby('rbs')

sns.violinplot(x='rbs', y='relative_intensity',
               data=data, linewidth=0.75, alpha=0.5)
