import glob
import skimage.io
import numpy as np
import skimage.io
import warnings

# Define the experiment paramters
DATE = 20171206
BASENAME = 'sfGFP_10ngmL'
missing_channel = 2
num_timepoints = 37
num_positions = 16
# ---------------------------------------------------------------------------
# Shouldn't need to change anything below here.
# ---------------------------------------------------------------------------
data_dir = '../../../data/images/{0}_{1}_dilution/growth'.format(
    DATE, BASENAME)
# Loop through each time point.
for i in range(num_positions):
    for j in range(num_timepoints):
        im_files = glob.glob(
            '{0}/*t{1:05d}*xy{2:03d}*c*.tif'.format(data_dir, j, i))
        chan_files = glob.glob(
            '{0}/*t{1:05d}*xy{2:03d}*c{3}.tif'.format(data_dir, j, i, missing_channel))

        if len(chan_files) == 0:
            # Get the new blank image name.
            ex_im = skimage.io.imread(im_files[0])
            blank_im = np.zeros_like(ex_im)
            new_name = '{0}{1}.tif'.format(im_files[0][:-5], missing_channel)

            # Ignore the "low contrast" warning.
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                skimage.io.imsave(new_name, blank_im)
