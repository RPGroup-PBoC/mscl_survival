import numpy as np
import matplotlib.pyplot as plt
import bokeh.io
import mscl.plotting
import mscl.stats
import mscl.process
import pandas as pd
import glob
import os
import matlab.engine as m
import imp
imp.reload(mscl.stats)
eng = m.start_matlab()
colors = mscl.plotting.set_plotting_style()
# Set the experimental parameters.
DATE = 20171129
BASE_STRAIN = 'sfGFP_10ngmL'
ATC_CONC = 10  # in ng/mL

# %%

# Define the directory containing the autofluorescence samples.
root_dir = '../../../data/images/{0}_{1}_dilution/'.format(DATE, BASE_STRAIN)

auto_dir = '{0}auto/xy*'.format(root_dir)
auto_positions = glob.glob(auto_dir)

# Loop through each position and parse the cell files.
auto_dfs = []
for i, pos in enumerate(auto_positions):
    # Get the positition number.
    num = int(pos.split('xy')[-1])

    # Glob the cell files and generkate the Datarame.
    files = glob.glob('{0}/cell/*.mat'.format(pos))
    _df = mscl.process.parse_cell_files(files, eng, add_props={'position': int(num)},
                                        excluded_props=['daughter_1_ID', 'daughter_2_ID',
                                                        'death_fluo', 'motherID',
                                                        'sisterID', 'divide', 'birth', 'death'])
    auto_dfs.append(_df)

# Assemble all of the positions into a single DataFrame
auto_df = pd.concat(auto_dfs, axis=0)
# Compute the mean autofluorescence value.
mean_auto = auto_df['birth_fluo'].mean()
print('The mean autofluorescence is {0:.2f} a.u. per cell.'.format(mean_auto))

# Plot and save the autofluorescence distribution.
x, y = mscl.stats.ecdf(auto_df['birth_fluo'])
fig, ax = plt.subplots(1, 1)
ax.set_xlabel('integrated autofluorescence [a.u.]')
ax.set_ylabel('ECDF')
_ = ax.step(x, y, color='dodgerblue')
_ = ax.vlines(mean_auto, 0, 1.0, color='tomato', lw=3, alpha=0.5)
plt.tight_layout()
plt.savefig('output/{0}_{1}ngmL_autofluorescence_distribution.png'.format(
    DATE, ATC_CONC), bbox_inches='tight')


# %% Process the growth data.
# Assemble the entire DataFrame.
data_dir = '{0}growth/xy*'.format(root_dir)
data_dir
# Grab all of the positions.
positions = glob.glob(data_dir)

# Loop through each position and generate the DataFrame.
dfs = []
for i, pos in enumerate(positions):
    print(i)
    # Figure out the position number.
    num = int(pos.split('xy')[-1])

    # Grab all of the files.
    files = glob.glob('{0}/cell/*.mat'.format(pos))

    # Generate the DataFrame.
    _df = mscl.process.parse_cell_files(
        files, eng, add_props={'position': num})
    dfs.append(_df)
dilution_df = pd.concat(dfs, axis=0)

#%% Plot and save the initial intensity disribution.
founders = dilution_df.loc[dilution_df['birth'] == 1]['birth_fluo']
x, y = mscl.stats.ecdf(founders)
fig, ax = plt.subplots(1, 1)
ax.set_xlabel('integrated dilution strain fluorescence [a.u.]')
ax.set_ylabel('ECDF')
_ = ax.step(x, y, color='dodgerblue')
plt.tight_layout()
plt.savefig('output/{0}_{1}ngmL_dilution_strain_distribution.png'.format(
    DATE, ATC_CONC), bbox_inches='tight')


# %% Check for conservation of fluorescence.
final_frame = dilution_df['death'].max()
# Restrict data frame to only those that divided or died on last frame.
producers = dilution_df[(dilution_df['divide'] == 1) | (
    dilution_df['death'] == final_frame)].copy()

# Remove anomolous zeros.
measured = producers.loc[(producers['birth_fluo'] > 0)
                         & (producers['death_fluo'] > 0)]
measured['birth_fluo'] -= mean_auto
measured['death_fluo'] -= mean_auto

grouped = measured.groupby(['position', 'motherID'])

# Iterate through and compute the total fluorescence and mother fluorescence.
summed = []
mother = []
for g, d in grouped:
    mother_fluo = (measured[(measured['ID'] == g[1]) & (
        measured['position'] == g[0])]['death_fluo']).sum()
    if mother_fluo != 0:
        mother.append(mother_fluo)
        summed.append(np.sum(d['birth_fluo']))

# Plot the conservation
# compute the expected linear trend.
theo = np.linspace(np.min(mother), np.max(mother), 300)
fig, ax = plt.subplots(1, 1)
ax.set_xlabel('mother fluorescence [a.u.]')
ax.set_ylabel('summed daughter fluorescence [a.u.]')
ax.plot(mother, summed, 'o', color='slategray', alpha=0.5, label='data', ms=1)
ax.plot(theo, theo, color='dodgerblue', lw=2, label='conserved fluorescence')
plt.legend()
plt.tight_layout()
plt.savefig('output/{0}_{1}_fluorescence_conservation.png'.format(DATE,
                                                                  BASE_STRAIN), bbox_inches='tight')
