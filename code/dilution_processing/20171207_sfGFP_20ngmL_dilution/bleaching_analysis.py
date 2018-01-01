import numpy as np
import matplotlib.pyplot as plt
import mscl.plotting
import mscl.image
import mscl.stats
import mscl.mcmc
import pandas as pd
import glob
import skimage.io
import theano.tensor as tt
from tqdm import tqdm
import pymc3 as pm
colors = mscl.plotting.set_plotting_style()

# Define the experimental parameters.
DATE = 20171207
BASE_STRAIN = 'sfGFP_20ngmL'
SAMPLE_NAMES = ['auto', 'delta']
root_dir = '../../../data/images/{0}_{1}_dilution/'.format(
    DATE, BASE_STRAIN)
IP_DIST = 0.16
EXPOSURE_MS = 100  # Exposure in ms.

# %% Load the flatfield correction.
slide_files = glob.glob('{0}/*fluorescent_slide*/*/*.tif'.format(root_dir))
slide_ims = skimage.io.ImageCollection(slide_files)
mean_field = np.mean(slide_ims, axis=0)


# %% Process the bleaching files.
col_names = ['strain', 'cell_id', 'area', 'total_intensity',
             'elapsed_time_s', 'mean_bg', 'date']
bleaching_df = pd.DataFrame([], columns=col_names)
print('Beginning processing...')
for i, nom in enumerate(tqdm(SAMPLE_NAMES)):
    files = glob.glob(
        '{0}photobleaching/*{1}_bleaching*'.format(root_dir, nom))
    cell_counter = 0
    for j, pos in enumerate(tqdm(files, desc='Processing positions')):
        # Load the images
        file = glob.glob('{0}/*.tif'.format(pos))
        ims = skimage.io.ImageCollection(file)

        # Split in to phase and intensity.
        bf_im = ims[0][:, 0, :, :][0]
        fluo_ims = ims[0][:, 1, :, :]

        # Segment the image and generate the inverse mask.
        mask = mscl.image.threshold_phase(bf_im)
        inv_mask = (mask < 1)

        # Flatten the fluorescence images.
        fluo_flat = [mscl.image.generate_flatfield(
            im, mean_field) for im in fluo_ims]
        for k, im in enumerate(fluo_flat):
            # Compute the mean background of the position.
            mean_bg = np.mean(im[inv_mask])
            props = skimage.measure.regionprops(mask, im)
            cell_no = 0
            for z, p in enumerate(props):
                area = p.area * IP_DIST**2
                total_intensity = p.mean_intensity * area
                label = z + cell_counter
                elapsed_time = k * EXPOSURE_MS / 1E3
                cell_dict = dict(strain=nom, cell_id=label,
                                 total_intensity=total_intensity - area * mean_bg,
                                 elapsed_time_s=elapsed_time,
                                 mean_bg=mean_bg, date=DATE, area=area)
                bleaching_df = bleaching_df.append(
                    cell_dict, ignore_index=True)
        cell_counter += np.max(mask)
print('...bleaching files processed, begining analysis...')

# %% Rescale data
grouped = bleaching_df.groupby(['strain', 'cell_id'])
dfs = []
for g, d in grouped:
    i0 = d[d['elapsed_time_s'] == 0.0]['total_intensity'].values[0]
    d.loc[:, 'rescaled_intensity'] = d.loc[:, 'total_intensity'] / i0
    dfs.append(d)
rescaled_df = pd.concat(dfs, axis=0)
rescaled_df.to_csv('output/{0}_sfGFP_bleaching.csv'.format(DATE))

# %% Compute the mean traces for each.
grouped = rescaled_df.groupby(['elapsed_time_s'])
mean_df = pd.DataFrame([], columns=['mean_auto', 'mean_delta',
                                    'elapsed_time', 'sub_mean'])

for g, d in grouped:
    mean_auto = d[d['strain'] == 'auto']['total_intensity'].mean()
    mean_delta = d[d['strain'] == 'delta']['total_intensity'].mean()
    sub = mean_delta - mean_auto
    cell_dict = dict(elapsed_time=g, mean_auto=mean_auto,
                     mean_delta=mean_delta, sub_mean=sub)
    mean_df = mean_df.append(cell_dict, ignore_index=True)

mean_df['sub_rescaled'] = mean_df['sub_mean'] / mean_df.iloc[0]['sub_mean']

# %% Perform inference of autofluorescence bleaching to set informative priors.
with pm.Model() as bleach_model:
    # Set uninformative priors.
    _tau_1 = pm.HalfNormal('tau_1', sd=100)
    _tau_2 = pm.HalfNormal('tau_2', sd=100)
    _beta_1 = pm.Uniform('beta_1', lower=0, upper=1)
    _beta_2 = pm.Uniform('beta_2', lower=0, upper=1)
    _bg = pm.Uniform('bg', lower=0, upper=1)
    sigma = pm.HalfNormal('sigma', sd=1)

    # Set the observables.
    time = mean_df['elapsed_time'].values
    obs = mean_df['sub_rescaled'].values

    # Compute the expected value.
    I_theo = _bg + _beta_1 * tt.exp(-time / _tau_1) +\
        _beta_2 * tt.exp(-time / _tau_2)

    # Set the likelihood and sample.
    like = pm.Normal('like', mu=I_theo, sd=sigma, observed=obs)
    bleach_trace = pm.sample(draws=5000, tune=1000)

# %%
mcmc_df = mscl.mcmc.trace_to_dataframe(bleach_trace, bleach_model)
fit_stats = mscl.mcmc.compute_statistics(mcmc_df)
fit_stats.to_csv('output/{0}_{1}_bleaching_constants.csv'.format(DATE, BASE_STRAIN),
                 index=False)
tau_1 = fit_stats[fit_stats['parameter'] == 'tau_1']['mode'].values[0]
tau_2 = fit_stats[fit_stats['parameter'] == 'tau_2']['mode'].values[0]
beta_1 = fit_stats[fit_stats['parameter'] == 'beta_1']['mode'].values[0]
beta_2 = fit_stats[fit_stats['parameter'] == 'beta_2']['mode'].values[0]
bg = fit_stats[fit_stats['parameter'] == 'bg']['mode'].values[0]

# Plot the fit.
time_range = np.linspace(0, 100, 1000)
fit = bg + beta_1 * np.exp(-time_range / tau_1) + \
    beta_2 * np.exp(-time_range / tau_2)
fig, ax = plt.subplots(1, 1)
ax.set_xlabel('time [s]')
ax.set_ylabel('rescaled intensity')
_ = ax.plot(mean_df['elapsed_time'], mean_df['sub_rescaled'], '.', color='slategray',
            label='auto subtracted mean')
_ = ax.plot(time_range, fit, color='dodgerblue', label='biexponential fit')
plt.legend()
plt.tight_layout()
plt.savefig('output/{0}_sfGFP_bleaching_fit.png'.format(DATE))
