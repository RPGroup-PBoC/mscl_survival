"""
Rilename: fig2.py
Author: Griffin Chure
Creation Date: 20180119
Purpose: To generate comparison images of the expression level of RBS mutants.
"""
import numpy as np
import matplotlib.pyplot as plt
import mscl.plotting
import mscl.image
import skimage.io
import glob
colors = mscl.plotting.set_plotting_style()


# Define the data directory.
data_dir = '../data/images/'

# Load the fluorescent slide image for ff correction. Note this is a different day.
field_files = glob.glob(
    '{0}20170922_sd0*calibration/Pos*/*GFP*.tif'.format(data_dir))
field_ims = skimage.io.ImageCollection(field_files)

# Create the mean projection.
field_im = mscl.image.average_stack(field_ims)

# Handpick some images for comparison.
SD4 = skimage.io.imread(
    '{0}20170714_sd4_pre_5ulmin_0.018hz_shock/SD4_5ulmin_pre_2_MMStack_Pos5.ome.tif'.format(data_dir))[1, :, :]
SD2 = skimage.io.imread(
    '{0}/20170517_sd2_pre_1200ulmin_2.17hz_shock/Pos1/img_000000000_GFP_000.tif'.format(data_dir))
SD2_12 = skimage.io.imread(
    '{0}/20170519_12sd2_pre_5ulmin_0.50hz_shock/Pos16/img_000000000_GFP_000.tif'.format(data_dir))
SD1 = skimage.io.imread(
    '{0}/20170525_sd1_pre_5ulmin_0.01hz_shock/Pos3/img_000000000_GFP_000.tif'.format(
        data_dir))

# Manually define the exposure times.
EXP_TIMES = {'SD4': 100, 'SD2': 50, '12SD2': 50, 'SD1': 25}

ims = {'SD4': SD4, 'SD2': SD2, '12SD2': SD2_12, 'SD1': SD1}
ff_ims = {}
for i in ims:
    im = ims[i]
    ff = mscl.image.generate_flatfield(im, field_im)
    ff_ims[i] = ff

# Correct the expsure for all images.
for i in ff_ims:
    ff_ims[i] = ff_ims[i] * EXP_TIMES['SD4'] / EXP_TIMES[i]

# Plot all four images together.
fig, ax = plt.subplots(1, 4)
min_int = np.min(ff_ims)
max_int = np.max(ff_ims)
for i, a in zip(ff_ims.keys(), ax):
    _ = a.imshow(ff_ims[i], cmap='viridis', vmin=20, vmax=2E4)
