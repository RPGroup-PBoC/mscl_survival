import numpy as np
import glob
import pandas as pd
import matlab.engine as matlab
try:
    eng = matlab.connect_matlab()
except:
    eng = matlab.start_matlab()

# Define the experiment parameters
DATE = 20171012
BASENAME = 'sfGFP_6ngmL'
STRAIN = 'HG105 galK::sfGFP'
ATC_CONC = 6

# Define the data directory for segmentation.
data_dir = '../../../data/images/{0}_{1}_dilution/'.format(DATE, BASENAME)

# %% Segment via MatLab SuperSegger.
status = eng.dilutionSegmentation(data_dir)
print(status)


# %% Process the cell files.
data_dir = '/Users/gchure/Desktop/seg_test'
files = glob.glob(data_dir + '*.mat')


vals = ['position', 'date', 'strain', 'atc_ngmL', 'death_frame', 'birth_frame',
        'division', 'ID', 'mother_ID', 'sister_ID', 'daughter_1_ID',
        'daughter_2_ID', 'birth_fluo',
        'death_fluo']

df = pd.DataFrame([], columns=vals)
for f in files:
    # Load the matrix using Matlab for "easy" parsing.
    eng.workspace['f'] = f
    mat = eng.eval("load(f)")

    # Extract ubiquitous properties.
    cell_dict = {'date': DATE, 'strain': STRAIN, 'atc_ngmL': ATC_CONC,
                 'death_frame': mat['death'], 'birth_frame': mat['birth'],
                 'division': bool(mat['divide']), 'ID': mat['ID'],
                 'mother_ID': mat['motherID'], 'sister_ID': mat['sisterID']}

    # Death with daughter ID issue.
    daughters = mat['daughterID']
    if len(daughters) is 0:
        cell_dict['daughter_1_ID'] = None
        cell_dict['daughter_2_ID'] = None
    else:
        cell_dict['daughter_1_ID'] = daughters[0][0]
        cell_dict['daughter_2_ID'] = daughters[0][1]

    # Parse the Cell Array for fluorescence information.
    cell_dict['birth_fluo'] = mat['CellA'][0]['fl1']['sum']
    cell_dict['death_fluo'] = mat['CellA'][-1]['fl1']['sum']

    # Append to the DataFrame
    df = df.append(cell_dict, ignore_index=True)
