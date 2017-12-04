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
eng = m.start_matlab()
colors = mscl.plotting.set_plotting_style()

# Set the experimental parameters.
DATE = 20171130
BASE_STRAIN = 'sfGFP_10ngmL'
EXPOSURE_MS = 100
ATC_CONC = 10


# ---------------------------------------------------------------------------
# PROCESS CELL OUTPUT FROM SUPERSEGGER
# -----------------------------------------------------------------------------
# %%  Generate dataframes of the cells
# Define the directory containing the autofluorescence samples.
root_dir = '../../../data/images/{0}_{1}_dilution/'.format(DATE, BASE_STRAIN)

auto_dir = '{0}autofluorescence/xy*'.format(root_dir)
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

# %% Process the growth data.
# Assemble the entire DataFrame.
data_dir = '{0}growth/xy*'.format(root_dir)

# Grab all of the positions.
positions = glob.glob(data_dir)

# Loop through each position and generate the DataFrame.
dfs = []
for i, pos in enumerate(positions):

    # Figure out the position number.
    num = int(pos.split('xy')[-1])

    # Grab all of the files.
    files = glob.glob('{0}/cell/*.mat'.format(pos))

    # Generate the DataFrame.
    _df = mscl.process.parse_cell_files(
        files, eng, add_props={'position': num})
    dfs.append(_df)
dilution_df = pd.concat(dfs, axis=0)


# ---------------------------------------------------------------------------
# CORRECT FOR PHOTOBLEACHING
# -----------------------------------------------------------------------------
#%%
bleach = pd.read_csv(
    'output/{0}_{1}_bleaching_constants.csv'.format(DATE, BASE_STRAIN))

# Extract the parameters.
bg = bleach[bleach['parameter'] == 'bg']['mode'].values[0]
tau_1 = bleach[bleach['parameter'] == 'tau_1']['mode'].values[0]
tau_2 = bleach[bleach['parameter'] == 'tau_2']['mode'].values[0]
beta_1 = bleach[bleach['parameter'] == 'beta_1']['mode'].values[0]
beta_2 = bleach[bleach['parameter'] == 'beta_2']['mode'].values[0]

# Correct the death fluorescence based on the number of exposures.
expo = dilution_df['num_exposures'] * EXPOSURE_MS / 1e3

bleaching_fraction = bg + beta_1 * \
    np.exp(-expo / tau_1) + beta_2 * np.exp(-expo / tau_2)
dilution_df.loc[:, 'corrected_fluo'] = (
    dilution_df.loc[:, 'death_fluo'] - mean_auto) / bleaching_fraction

# Correct for the mean autofluorescence.

# Save the dataframe.
dilution_df.to_csv('output/{0}_{1}_lineages.csv'.format(DATE, BASE_STRAIN))

#%% Plot and save the initial intensity disribution.
auto_x, auto_y = mscl.stats.ecdf(auto_df['birth_fluo'])
founders = dilution_df.loc[dilution_df['birth']
                           == 1]['corrected_fluo'] + mean_auto
dil_x, dil_y = mscl.stats.ecdf(founders)
fig, ax = plt.subplots(1, 1)
ax.set_xlabel(' intensity [a.u.]')
ax.set_ylabel('ECDF')
_ = ax.step(auto_x, auto_y, color='dodgerblue', label='autofluorescence')
_ = ax.vlines(mean_auto, 0, 1.0, color='tomato', lw=3, alpha=0.5,
              label='mean autofluorescence value')
_ = ax.step(dil_x, dil_y, color='tomato', label='dilution distribution')
plt.legend()
plt.tight_layout()
plt.savefig('output/{0}_{1}_dilution_strain_distribution.png'.format(
    DATE, BASE_STRAIN), bbox_inches='tight')


# -----------------------------------------------------------------------------
# CHECK CONSERVATION OF FLUORESCENCE
# -----------------------------------------------------------------------------
# %%
final_frame = dilution_df['death'].max()

# Restrict data frame to only those that divided or died on last frame.
producers = dilution_df[(dilution_df['divide'] == 1) | (
    dilution_df['death'] == final_frame)].copy()

# Remove anomolous zeros.
measured = producers.loc[(producers['birth_fluo'] > 0)
                         & (producers['death_fluo'] > 0)]

grouped = measured.groupby(['position', 'motherID'])

# Iterate through and compute the total fluorescence and mother fluorescence.
summed = []
mother = []
for g, d in grouped:
    mother_fluo = (measured[(measured['ID'] == g[1]) & (
        measured['position'] == g[0])]['corrected_fluo']).sum()
    if mother_fluo != 0:
        mother.append(mother_fluo)
        summed.append(np.sum(d['corrected_fluo']))

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


# -----------------------------------------------------------------------------
# ESTIMATE CALIBRATION FACTOR
# -----------------------------------------------------------------------------
# %%

# Remove cells with zero fluorescence intensity.
measured = producers.loc[producers['corrected_fluo'] > 0]

# Group sisters
children = measured.groupby(['position', 'motherID'])

# Instantiate the storage dataframe.
sister_df = pd.DataFrame(
    [], columns=['date', 'position', 'I1', 'I2', 'summed_int', 'sq_diff', 'cal_factor',
                 'cal_sd'])
for g, d in children:
    if len(d) == 2:
        I1, I2 = d['corrected_fluo'].values
        summed = I1 + I2
        sq_diff = (I1 - I2)**2

        # Assemble the dict and append to the dataframe.
        pair_dict = dict(date=DATE, position=g[0], I1=I1, I2=I2, summed_int=summed,
                         sq_diff=sq_diff, cal_factor=None)
        sister_df = sister_df.append(pair_dict, ignore_index=True)

# Estimate the calibration factor and sd.
I1, I2 = sister_df[['I1', 'I2']].values.T
popt = mscl.stats.estimate_calibration_factor(I1, I2)
sister_df.loc[:, ['cal_factor', 'cal_sd']] = popt
sister_df.to_csv(
    'output/{0}_{1}_calibration_factor.csv'.format(DATE, BASE_STRAIN), index=False)

# Generate the fit.
min_sum = np.round(np.log10(sister_df['summed_int'].min()))
max_sum = np.round(np.log10(sister_df['summed_int'].max()))
Itot_range = np.logspace(min_sum - 0.5, max_sum + 0.5, 500)
fit = popt[0] * Itot_range
extrema = [(popt[0] + popt[1]) * Itot_range, (popt[0] - popt[1]) * Itot_range]

# Bin the data for a sanity check.
bin_size = 50
sorted_vals = sister_df.sort_values(by=['summed_int'])
mean_sum = []
mean_sqdiff = []
bins = np.arange(0, len(sorted_vals) + bin_size, bin_size)
for i in range(1, len(bins)):
    d = sorted_vals.iloc[bins[i - 1]: bins[i]]
    mean_sum.append(d['summed_int'].mean())
    mean_sqdiff.append(d['sq_diff'].mean())

# Set the figure canvas.
fig, ax = plt.subplots(1, 1, figsize=(6, 4))
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel("$I_1 + I_2$ [a.u.]")
ax.set_ylabel('$(I_1 - I_2)^2$ [a.u.]')

# Plot the raw data
_ = ax.plot(sister_df['summed_int'], sister_df['sq_diff'], '.', color='slategray',
            alpha=0.5, label='raw data')

# Plot the binned data.
_ = ax.plot(mean_sum, mean_sqdiff, '.', color='tomato',
            label='binned data\n ($N$={0})'.format(bin_size))

# Plot the fit
_ = ax.plot(Itot_range, fit, color='dodgerblue', label='fit', lw=1)
_ = ax.fill_between(Itot_range, extrema[0], extrema[1], color='dodgerblue',
                    alpha=0.5, label='__nolegend__')
_ = ax.legend()
ax.set_title('α = {0} ± {1} a.u. / mol.'.format(int(popt[0]), int(popt[1])))
plt.tight_layout()
plt.savefig('output/{0}_{1}_calibration_factor.png'.format(DATE, BASE_STRAIN),
            bbox_inches='tight')
