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

# %%
%matplotlib inline
fig, ax = plt.subplots(1, 1)
ax.loglog(dilution_df['summed_int'], dilution_df['square_difference'],
          'o', color='slategray', alpha=0.5)


#%%
# Try to link up conservation of fluorescence from birth to death.
# Take only one position.
xy01 = growth_df.loc[growth_df['position'] == '01']

# Look only at cells that divided or died on last frame.
producers = xy01.loc[(xy01['division'] == True) & (xy01['death_frame'] != 26)]

#  | ((xy01['death_frame']==26) & (xy01['birth_frame'] != 1))]
last_generation = xy01[(xy01['death_frame'] == 26) &
                       (xy01['birth_frame'] != 1)]

# Insert an empty tree number.
producers.loc[:, 'tree_number'] = np.nan
last_generation.loc[:, 'tree_number'] = np.nan


# Set the tree roots.
producers.loc[producers['mother_ID'] == 0, 'tree_number'] = producers['ID']

# %% it_counter = 0
it_counter = 0
orphans = producers.isnull().values.any()
while (orphans == True) & (it_counter < 20):
    # Grab all of the NaN and non NaN values.
    unassigned = producers[producers.isnull().any(axis=1)]

    num_before = len(unassigned)
    assigned = producers[~producers.isnull().any(axis=1)]
    grouped = assigned.groupby('tree_number')
    for g, d in grouped:
        # Get a list of the daughter IDs.
        daughters = d[['daughter_1_ID', 'daughter_2_ID']].values.flatten()

        # Change the tree number in the producers. for those indices.
        for progeny in daughters:
            producers.loc[producers['ID'] == progeny, 'tree_number'] = g
    orphans = producers.isnull().values.any()
    orphans
    print(it_counter)
    it_counter += 1


# Link the final frame.
grouped == producers.groupby('tree_number')
for g, d in grouped:
    daughters = d[['daughter_1_ID', 'daughter_2_ID']].values.flatten()
    for progeny in daughters:
        last_generation.loc[last_generation['ID']
                            == progeny, 'tree_number'] = g

lineages = pd.concat([producers, last_generation], axis=0)

print('finished!')
# %% I'm skeptical this is working with just a single iteration.
# Ignore the ends for now.
complete_lineages = lineages.dropna()
complete_lineages.head()
# Group by the lineage information.
family_tree = complete_lineages.groupby('tree_number')
for g, d in family_tree:
    d

#
