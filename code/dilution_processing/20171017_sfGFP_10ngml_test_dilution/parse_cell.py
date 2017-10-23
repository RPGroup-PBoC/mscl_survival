import numpy as np
import glob
import pandas as pd
import mscl_utils as mscl
import matplotlib.pyplot as plt
import matlab.engine as matlab
import scipy.io
try:
    eng = matlab.connect_matlab()
except:
    eng = matlab.start_matlab()
colors = mscl.set_plotting_style()
# Define the experiment parameters
DATE = 20171017
BASENAME = 'sfGFP_10ngmL'
STRAIN = 'HG105 galK::sfGFP'
ATC_CONC = 10

# Define the data directory for segmentation.
data_dir = '../../../data/images/{0}_{1}_dilution/'.format(DATE, BASENAME)

# %%


def parse_cell_files(file_list):
    """
    Parses a directory of Cell mat files and returns a tidy Pandas DataFrame.
    """

    # Define the properties I care about.
    vals = ['position', 'date', 'strain', 'atc_ngmL', 'death_frame',
            'birth_frame', 'division', 'ID', 'mother_ID', 'sister_ID',
            'daughter_1_ID', 'daughter_2_ID', 'birth_fluo', 'death_fluo',
            'status']

    df = pd.DataFrame([], columns=vals)
    for i, f in enumerate(file_list):
        # Load the matrix through matlab for "easy" parsing.
        eng.workspace['f'] = f
        mat = eng.eval('load(f)')

        # Extract the standard properties
        cell_dict = {'position': p.split('xy')[-1], 'date': DATE, 'strain':
                     STRAIN, 'atc_ngmL': ATC_CONC, 'death_frame':
                     mat['death'], 'birth_frame':
                     mat['birth'], 'division': bool(mat['divide']), 'ID':
                     mat['ID'], 'mother_ID': mat['motherID'], 'sister_ID':
                     mat['sisterID'], 'status': mat['stat0']}

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
    return df


# %% Parse the auto fluorescence images.
positions = glob.glob('{0}autofluorescence/xy*'.format(data_dir))
auto_dfs = []
for i, p in enumerate(positions):
    files = glob.glob('{0}/cell/*.mat'.format(p))
    eng.workspace['f'] = files[0]
    mat = eng.eval('load(f)')
    _df = parse_cell_files(files)
    auto_dfs.append(_df)
auto_df = pd.concat(auto_dfs, ignore_index=True)

# Compute the mean autofluroescence value.
mean_auto = auto_df['birth_fluo'].mean()

# Parse the clist for the experiment directories.
positions = glob.glob('{0}growth/xy*'.format(data_dir))
growth_df = []
for i, p in enumerate(positions):
    files = glob.glob('{0}/cell/*.mat'.format(p))
    _df = parse_cell_files(files)
    growth_df.append(_df)
growth_df = pd.concat(growth_df, axis=0)

# %% Plot the fluctuations.
final_cells = growth_df[growth_df['death_fluo'] > 0]
grouped = final_cells.groupby(['mother_ID', 'position'])
dilution_df = pd.DataFrame([], columns=['square_difference', 'summed_int'])
for g, d in grouped:
    if len(d) == 2:
        sq_diff = np.diff(d['death_fluo'] - mean_auto)**2
        summed = (d['death_fluo'] - mean_auto).sum()
        cell_dict = {'square_difference': sq_diff,
                     'summed_int': summed}
        dilution_df = dilution_df.append(cell_dict, ignore_index=True)
# Save the dataframes.
auto_df.to_csv('../../../data/{0}_auto_supersegger_output.csv', index=False)
growth_df.to_csv('../../../data/{0}_supersegger_output.csv', index=False)
# %%
%matplotlib inline
fig, ax = plt.subplots(1, 1)
ax.loglog(dilution_df['summed_int'], dilution_df['square_difference'],
          'o', color='slategray', alpha=0.5)

# %%
# Try to link up conservation of fluorescence from birth to death.
# I only want to track the
xy01 = growth_df.loc[growth_df['position'] == '01']

# Include only the cells that divided or died on the final frame.
death_numbers = xy01['death_frame'].unique()
max_frame = np.max(death_numbers)
producers = xy01.loc[(xy01['division'] == True) &
                     (xy01['death_frame'] != max_frame)]

# Set the tree group number to NaN.
producers.loc[:, 'tree_number'] = np.nan
# Set the tree number to the mother ID.
producers.loc[producers['mother_ID'] == 0, 'tree_number'] = producers['ID']
producers[producers['mother_ID'] == 0]
# Set up a while loop that continues until there are no more nans.
orphans = producers.isnull().values.any()
iterations = 0
while orphans == True & iterations < 20:
    # Split the producerts into assigned and unassigned.
    assigned = producers[~producers['tree_number'].isnull()]
    unassigned = producers[producers['tree_number'].isnull()]
    len(producers)
    len(assigned)
    len(unassigned)
    unassigned
