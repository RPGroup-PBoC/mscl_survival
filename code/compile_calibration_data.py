#!/usr/bin/env python
import pandas as pd
import glob


def main():
    # Scrape the directories.
    root_dir = './calibration_processing/*/output/*.csv'
    paths = glob.glob(root_dir)

    # Screen through the list of ignored datasets.
    with open('excluded_datasets.txt', 'r') as f:
        ignored_files = f.read().splitlines()
        approved_paths = []
        for i, f in enumerate(paths):
            filename = f.split('/')[-1]
            if filename not in ignored_files:
                approved_paths.append(f)

    # Concatenate all dataframes together.
    df = (pd.read_csv(f, comment='#') for f in approved_paths)
    df = pd.concat(df, ignore_index=True)
    # Make new rescaled intensity column based off of the longest exposure.
    df['rescaled_intensity'] = (df['intensity'] - df['mean_bg']) *\
        df['exposure_ms'].max() / df['exposure_ms']
    df.to_csv('../data/csv/compiled_calibration_data.csv', index=False)


if __name__ == '__main__':
    main()
    print('All calibration data compiled. Thank you come again.')
