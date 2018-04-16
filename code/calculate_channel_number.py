# %%
# -*- coding: utf-8 -*-
import numpy as np
import glob
import pandas as pd

# Get the files.
files = glob.glob('processing/*shock*/output/*.csv')

# Read and concatenate.
csv = pd.concat([pd.read_csv(f, comment='#')
                 for f in files], ignore_index=True)

# Subtract the background and rescale the intensity.
max_exposure = csv['exposure_ms'].max()
csv.loc[:, 'rescaled_intensity'] = (csv['intensity'] - csv['mean_bg']) * max_exposure / csv['exposure_ms']


# Load the calibration factor.
cal_factor = pd.read_csv('../data/csv/calibration_factor.csv')

# Rescale the intensity and compute the mean.
<<<<<<< HEAD
calibration = (cal_factor['au_per_channel'] * max_exposure / cal_factor['exposure_ms']).values[0]
calibration_err = (cal_factor['au_per_channel_err'] * max_exposure / cal_factor['exposure_ms']).values[0]

# Compute the reference cell area.
avg_area = cal_factor['mean_area'].values[0]
err_area = cal_factor['mean_area_err'].values[0]

# Compute the channels per cell.
csv.loc[:, 'effective_channels'] = csv['rescaled_intensity'] * avg_area / calibration
csv.loc[:, 'effective_channels_err'] = csv['effective_channels'] * ((err_area / avg_area) + (calibration_err/calibration))
=======
calibration = cal_factor['au_per_channel'] * max_exposure / cal_factor['exposure_ms']
print(cal_factor)
mean_cal = np.mean(calibration)
sem_cal = np.std(calibration) / np.sqrt(len(calibration))

# Compute the reference cell area.
avg_area = cal_factor['area'].mean()


# Compute the channels per cell.
csv.loc[:, 'effective_channels'] = csv['rescaled_intensity'] * avg_area / mean_cal
>>>>>>> 709a6946ca985c61e256071636766d53e0e2ee5d

csv.loc[csv['flow_rate'] < 1.0, 'shock_class'] = 'slow'
csv.loc[csv['flow_rate'] >= 1.0, 'shock_class'] = 'fast'

csv.to_csv('../data/csv/compiled_data.csv', index=False)
