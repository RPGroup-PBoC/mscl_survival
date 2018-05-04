# -*- coding: utf-8 -*-
# %%
import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
sys.path.insert(0, '../../')
import mscl.stats
import mscl.plotting
colors = mscl.plotting.set_plotting_style()

# Load the data
data = pd.read_csv('../../data/csv/mscl_survival_data.csv')

# Generate a histogram and ECDF
fig, ax = plt.subplots(1, 2, figsize=(6, 4))

ax[0].hist(data['flow_rate'], bins=10)