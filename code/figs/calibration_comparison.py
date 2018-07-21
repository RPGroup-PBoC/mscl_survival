# %%
# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import pystan
import glob
import matplotlib.pyplot as plt
import scipy.stats
import skimage.io
import skimage.morphology
import skimage.measure
import skimage.filters
sys.path.insert(0, '../../')
import mscl.image
import mscl.stats
import mscl.plotting
import mscl.mcmc
colors = mscl.plotting.set_plotting_style()

#%% Define the data directory.
data_dir = '../../data/images/20180718_mlg910_comparison'

# Load the slide images and generate the average.
slide_files = glob.glob('{}/slide_2/*/*.tif'.format(data_dir))
avg_field = np.mean(skimage.io.ImageCollection(slide_files), axis=0)

# Define the media sources.
media = ['LB_miller', 'M9_glucose']
media_key = {'LB_miller': 1, 'M9_glucose': 2}
strains = ['MLG910']

# Loop through each media type
medium_df = pd.DataFrame([])
for i, m in enumerate(media):
    for j, st in enumerate(strains):
        bf_files = np.sort(glob.glob('{}/{}_{}_1/*/*Brightfield*.tif'.format(data_dir, m, st)))
        gfp_files = np.sort(glob.glob('{}/{}_{}_1/*/*GFP*.tif'.format(data_dir, m, st)))
        bf = skimage.io.ImageCollection(bf_files)
        gfp = skimage.io.ImageCollection(gfp_files)

        # Loop through each image.
        for k, im in enumerate(bf): 
            # Segment
            seg = mscl.image.contour_seg(im, level=0.2, min_int=0.4)

            # Flatten the fluorescence image.
            flat = mscl.image.generate_flatfield(gfp[k], avg_field)

            # Compute the mean background.
            mean_bg = mscl.image.compute_mean_bg(im, flat)

            # Extract properties
            props = skimage.measure.regionprops(seg, flat)
            for z, p in enumerate(props):
                area = p.area * 0.160**2
                mean_int = (p.mean_intensity - mean_bg) / 0.160**2
                cell_dict = dict(area=area, mean_int=mean_int, medium=m, strain=st,
                                 media_id=media_key[m])
                medium_df = medium_df.append(cell_dict, ignore_index=True)


# %% Load the stan model.
model = pystan.StanModel('../stan/calibration_factor.stan')

# %%

# %%
# Set up the data dictionary and sample.
medium_df = medium_df[medium_df['strain']=='MLG910'].copy()
medium_df[medium_df['medium'] == 'LB_miller']['area'] *= 0.5
data_dict = dict(J=2, N=len(medium_df), media=medium_df['media_id'].values.astype(int),
                mean_intensity=medium_df['mean_int'].values * 4, channel_mu=np.array([340, 466]),
                channel_sig=np.array([68, 64]), area=medium_df['area'].values)

data_dict
# %%
# Sample the model.
print('Sampling....')
samples = model.sampling(data_dict, iter=100000, chains=4, thin=50)
print('finished!')

# %%
samples
# %%
samp_df = mscl.mcmc.chains_to_dataframe(samples)
samp_df.to_csv('../../data/csv/calibration_factor_media_comparison.csv', index=False)
# %%
plt.hist(samp_df[samp_df['alpha_mu__0'] > 1000]['alpha_mu__0'], color=colors['blue'], alpha=0.5, bins=100,
    label='LB Miller')
plt.hist(samp_df['alpha_mu__1'], color=colors['red'], alpha=0.5, bins=100, label='M9 + 0.5% glucose')
plt.legend(fontsize=8)
plt.set



# %% Process the LB 500 data.
lb500_df = pd.DataFrame([])
bf_files = np.sort(glob.glob('{}/LB_500_MLG910_1/*/*Brightfield*.tif'.format(data_dir)))
gfp_files = np.sort(glob.glob('{}/LB_500_MLG910_1/*/*GFP*.tif'.format(data_dir, m, st)))
bf = skimage.io.ImageCollection(bf_files)
gfp = skimage.io.ImageCollection(gfp_files)

# Loop through each image.
for k, im in enumerate(bf): 
    # Segment
    seg = mscl.image.contour_seg(im, level=0.2, min_int=0.4)

    # Flatten the fluorescence image.
    flat = mscl.image.generate_flatfield(gfp[k], avg_field)

    # Compute the mean background.
    mean_bg = mscl.image.compute_mean_bg(im, flat)

    # Extract properties
    props = skimage.measure.regionprops(seg, flat)
    for z, p in enumerate(props):
        area = p.area * 0.160**2
        mean_int = (p.mean_intensity - mean_bg) / 0.160**2
        cell_dict = dict(area=area, mean_int=mean_int, medium=m, strain=st,
                                 media_id='LB_500mM')
        lb500_df = lb500_df.append(cell_dict, ignore_index=True)

# %% Instantiate the figure
fig, ax = plt.subplots(1, 2, figsize=(5, 3))
# Format the axes.
for a in ax:
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)

# Add labels.
ax[0].set_xlabel('areal intensity [a.u. / µm$^2$]', fontsize=8)
ax[0].set_ylabel('cumulative distribution', fontsize=8)
ax[1].set_ylabel('estimated density', fontsize=8)
ax[1].set_xlabel('calibration factor [a.u. / channel]', fontsize=8)
ax[0].text(-0.3, 1.0, '(A)', fontsize=8, transform=ax[0].transAxes)
ax[1].text(-0.32, 1.0, '(B)', fontsize=8, transform=ax[1].transAxes)
# Compute the ecdfs
glucose_x, glucose_y = mscl.stats.ecdf(medium_df[medium_df['medium']=='M9_glucose']['mean_int'])
miller_x, miller_y = mscl.stats.ecdf(medium_df[medium_df['medium']=='LB_miller']['mean_int'])
lb500_x, lb500_y = mscl.stats.ecdf(lb500_df['mean_int'])

# Plot the distributions.
ax[0].step(glucose_x, glucose_y, lw=2, color=colors['red'], label='M9 +\n0.5% glucose')
ax[0].step(miller_x, miller_y, lw=2, color=colors['blue'], label='LB Miller')
ax[0].step(lb500_x, lb500_y, lw=2, color=colors['green'], label='LB +\n500mM NaCl')
ax[0].legend(loc='lower right', fontsize=8, handlelength=1)
# Plot the calibration factor KDEs.
alpha_range = np.linspace(1000, 8000, 500)
miller_alpha_kernel = scipy.stats.gaussian_kde(samp_df['alpha_mu__0'])
miller_kde = miller_alpha_kernel(alpha_range)
miller_kde *= np.sum(miller_kde)**-1
glucose_alpha_kernel = scipy.stats.gaussian_kde(samp_df['alpha_mu__1'])
glucose_kde = glucose_alpha_kernel(alpha_range)
glucose_kde *= np.sum(glucose_kde)**-1
ax[1].plot(alpha_range, glucose_kde, color=colors['red'], label='M9 +\n0.5% glucose')
ax[1].fill_between(alpha_range, glucose_kde, color=colors['red'], label='__nolegend__', alpha=0.5)
ax[1].plot(alpha_range, miller_kde, color=colors['blue'], label='LB Miller')
ax[1].fill_between(alpha_range, miller_kde, color=colors['blue'], label='__nolegend__', alpha=0.5)
ax[1].legend(fontsize=8, handlelength=1)
plt.tight_layout()
plt.savefig('../../figs/figRX_calibration_factor.pdf')
plt.savefig('../../figs/figRX_calibration_factor.png', dpi=300)

