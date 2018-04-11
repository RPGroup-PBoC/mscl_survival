# -*- coding: utf-8 -*-
import numpy as np
import glob
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import sys
sys.path.insert(0, '../../../')
import mscl.plotting
colors = mscl.plotting.set_plotting_style()

DATE = 20180410
CHANNEL_NO = 340
RBS = 'mlg910'

# Load the data.
data = pd.read_csv(
    'output/{}_{}_intensity_calibration.csv'.format(DATE, RBS), comment='#')

# Group by the replicate number
grouped = data.groupby('replicate_number').mean()
mean_cell = grouped['area'] * (grouped['intensity'] - grouped['mean_bg']) / 2
cal_factor = mean_cell / CHANNEL_NO

# df = pd.DataFrame(grouped).reset_index()
# df.loc[:, 'corrected_intensity'] = mean_cell
# df.loc[:, 'au_per_channel'] = cal_factor
df = pd.DataFrame(np.array([grouped['area'].values / 2, grouped['exposure_ms'],
                            cal_factor]).T, columns=['area', 'exposure_ms', 'au_per_channel'])
df.to_csv('../../../data/csv/calibration_factor.csv', index=False)

# %% Compare with others.
old_data = pd.read_csv(
    '../20170721_mlg910_calibration/output/20170721_mlg910_intensity_calibration.csv', comment='#')
old_data.loc[:, 'replicate_number'] = 6
merged = pd.concat([data, old_data], ignore_index=True)
max_exp = merged['exposure_ms'].max()
merged.loc[:, 'intensity'] = (merged['intensity'] - merged['mean_bg']) * max_exp / merged['exposure_ms']

fig, ax = plt.subplots(1, 1)
sns.boxplot(merged['replicate_number'], merged['intensity'], ax=ax, linewidth=1)

ax.set_ylabel('intensity')
ax.set_xticklabels(['pad 1', 'pad 2', 'pad 3', 'pad 4', 'pad 5', 'flow cell'])

plt.tight_layout()
plt.savefig('output/intensity_comparison.png', bbox_inches='tight')
