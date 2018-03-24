#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
This script loops through all processed data files and compiles individual
csv files into a single csv.
"""
import numpy as np
import glob
import pandas as pd
import os


def compile():
    files = glob.glob('./processing/*/output/*.csv')
    dfs = []
    for f in files:
        if 'shock' in f:
            date, strain, shock, _ = f.split(
                './processing/')[1].split('/')[0].split('_')
            shock = float(shock.split('hz')[0])
            id = 'shock'
        elif 'calibration' in f:
            date, strain, _ = f.split(
                './processing/')[1].split('/')[0].split('_')
            id = 'calibration'
            shock = np.nan
        if shock < 1.0:
            shock_class = 'slow'
        elif shock >= 1.0:
            shock_class = 'fast'
        else:
            shock_class = np.nan

        df = pd.read_csv(f, comment='#')
        df.insert(0, 'shock_class', shock_class)
        df.insert(0, 'experiment', id)
        dfs.append(df)

    compiled_data = pd.concat(dfs)

    # Correct for exposure time
    exposure_correction = compiled_data['exposure_ms'].max(
    ) / compiled_data['exposure_ms']
    compiled_data.loc[:, 'scaled_intensity'] = (compiled_data['intensity'] -
                                                compiled_data['mean_bg']) *\
        exposure_correction

    # Compute the calibration factor.
    mlg910 = compiled_data[compiled_data['rbs'] == 'mlg910']
    CHANNEL_NUM = 340
    CHANNEL_ERR = 68
    PERCENT_ERR = CHANNEL_NUM / CHANNEL_ERR
    AU_PER_CHAN = np.mean(
        mlg910['area'] * mlg910['scaled_intensity']) / CHANNEL_NUM
    AVG_AREA = mlg910['area'].mean()

    compiled_data.loc[:, 'effective_channels'] = (AVG_AREA *
                                                  compiled_data['scaled_intensity']) / AU_PER_CHAN

    compiled_data.loc[:, 'effective_channel_err'] = PERCENT_ERR * \
        compiled_data['effective_channels']
    compiled_data = compiled_data[compiled_data['effective_channels'] > 0]

    # Drop unnecessary columns.
    compiled_data.drop(['intensity', 'mean_bg', 'replicate_number'],
                       axis=1, inplace=True)
    compiled_data.to_csv('../data/csv/mscl_survival_data.csv', index=False)


if __name__ == '__main__':
    compile()
