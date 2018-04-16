# -*- coding: utf-8 -*-
#%%
import numpy as np
import glob
import skimage.io
import bokeh.io
import bokeh.plotting
import pandas as pd
import os
import tqdm
import sys
sys.path.insert(0, '../../../')
import skimage.io
import scipy.ndimage
import mscl.plotting
import mscl.image


# %% Define the relevant parameters
DATE = 20180410
RBS = 'mlg910'
REP = '0000'
IP_DIST = 0.16 # µm / pixel

# Load in the data.
data_dir = '../../../data/images/{}_{}_pre_{}_intensity_calibration'.format(DATE, RBS, REP)

# Load in image files. 
mlg910_files = np.sort(glob.glob(data_dir + '/*ome.tif'))
metadata = np.sort(glob.glob(data_dir + '/*.txt'))
fluo_slide_files = glob.glob(data_dir + '/fluorescent_slide/*.ome.tif')

# Generate the average fluorescent slide.
fluo_slide_ims = skimage.io.ImageCollection(fluo_slide_files, conserve_memory=False)

# Median filter the images.
selem = skimage.morphology.square(3)
filt_fluo_slide = [scipy.ndimage.median_filter(i, footprint=selem) for i in fluo_slide_ims]
avg_field =  np.mean(filt_fluo_slide, axis=0)

# Choose a random image to save a segmentaiton mask.
seg_choice = np.random.choice(mlg910_files)
# Loop through each image. 
dfs = []
for i, f in enumerate(mlg910_files):
    # Get the identifiers.
    _, replicate, _, _, pos = f.split('/')[-1].split('_')
    pos = int(pos.split('.')[0].split('Pos')[1])

    # Load the image.
    im = skimage.io.imread(f)
    phase_im = im[0]
    fluo_im = im[1]

    # Flatten the field
    im_flat = mscl.image.generate_flatfield(fluo_im, avg_field)

    # Segment the image.
    im_seg = mscl.image.contour_seg(phase_im, level=0.2, min_int=0.3, perim_bounds=(5, 1E9), ecc_bounds=(0, 1),
    area_bounds=(0, 1E3))

    if f == seg_choice:
        mscl.plotting.save_seg('output/{}_{}_{}_example_segmentation.png'.format(DATE, RBS, replicate),
        phase_im, im_seg)        

    if np.max(im_seg) != 0:
        # Compute the mean background fluorescence.
        mean_bg = mscl.image.compute_mean_bg(phase_im, im_flat, method='isodata')

        # Extract the various properties.
        props = skimage.measure.regionprops(im_seg, im_flat)
        area = np.array([p.area * IP_DIST**2 for p in props])
        eccentricity = np.array([p.eccentricity for p in props])
        mean_int =  np.array([p.mean_intensity / IP_DIST**2 for p in props])

        # Scrape the metadata
        exposure = mscl.image.scrape_metadata(metadata[i], return_date=False)
        exposure = exposure['GFP_exp_ms']


        # Assemble the datframe.
        df = pd.DataFrame(np.array([area, eccentricity, mean_int]).T,\
                         columns=['area', 'eccentricity', 'intensity'])
        df.loc[:, 'date'] = DATE
        df.loc[:, 'rbs'] = RBS
        df.loc[:, 'mean_bg'] = mean_bg / IP_DIST**2
        df.loc[:, 'exposure_ms'] = exposure
        df.loc[:, 'replicate_number'] = int(replicate)
        dfs.append(df)    

# Concatenate the data frame.
data  = pd.concat(dfs, ignore_index=True)

target = 'output/{}_{}_intensity_calibration.csv'.format(DATE, RBS)

data.to_csv(target, index=False)



