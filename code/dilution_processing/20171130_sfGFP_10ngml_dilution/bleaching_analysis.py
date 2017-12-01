import numpy as np
import matplotlib.pyplot as plt
import mscl.plotting
import mscl.image
import mscl.stats
import mscl.mcmc
import pandas as pd
import glob
import skimage.io
from tqdm import tqdm
import pymc3 as pm

# Define the experimental parameters.
DATE = 20171129
BASE_STRAIN = 'sfGFP_10ngmL'
SAMPLE_NAMES = ['autofluorescence', 'delta']
root_dir = '../../../data/images/{0}_{1}_dilution/'.format(
    DATE, BASE_STRAIN)
IP_DIST = 0.16
EXPOSURE_MS = 100  # Exposure in ms.

# %% Load the flatfield correction.
slide_files = glob.glob('{0}/*fluorescent_slide*/*/*.tif'.format(root_dir))
slide_ims = skimage.io.ImageCollection(slide_files)
mean_field = np.mean(slide_ims, axis=0)


# %% Process the bleaching files.
col_names = ['strain', 'cell_id',
             'total_intensity', 'elapsed_time_s', 'mean_bg', 'date']
bleaching_df = pd.DataFrame([], columns=col_names)
print('Beginning processing...')
for i, nom in enumerate(tqdm(SAMPLE_NAMES)):
    files = glob.glob(
        '{0}photobleaching/{1}_bleaching*'.format(root_dir, nom))
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
                                 total_intensity=total_intensity,
                                 elapsed_time_s=elapsed_time,
                                 mean_bg=mean_bg, date=DATE)
                bleaching_df = bleaching_df.append(
                    cell_dict, ignore_index=True)
        cell_counter += np.max(mask)
print('...bleaching files processed, begining analysis...')

# %%
bleaching_df
