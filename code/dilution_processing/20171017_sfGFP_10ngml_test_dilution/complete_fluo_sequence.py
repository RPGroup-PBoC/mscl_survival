import glob
import skimage.io
import numpy as np
import skimage.io
import warnings

# Define the experiment paramters
DATE = 20171017
BASENAME = 'sfGFP_10ngmL'
missing_channel = 2
num_timepoints = 26
num_positions = 14
fluorescent_slide = True
num_auto_positions = 10
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


#%% Generate an average of the slide images.
slide_dir = data_dir.split('growth')[0]
dest_dir = [data_dir, '{0}autofluorescence'.format(slide_dir)]


if fluorescent_slide is True:
    # This will generate an average image of the illumination.
    slide_files = glob.glob(
        '{0}{1}_fluorescent_slide/Pos*/*.tif'.format(slide_dir, DATE))

    # Load in all files.
    im = skimage.io.ImageCollection(slide_files)
    bit_depth = im[0].dtype

    # Compute the mean image.
    mean_im = np.zeros_like(im[0])
    for i in im:
        mean_im += i
    mean_im = mean_im / len(im)
    mean_im = mean_im.astype(bit_depth)

    # Save the slide file repeatedly as the second fluorescent channel.
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        for d in dest_dir:
            if 'autofluorescence' in d:
                d
                bn = BASENAME + '_autofluorescence'
                num_pos = num_auto_positions
                num_time = 1
            else:
                bn = BASENAME
                num_pos = num_positions
                num_time = num_timepoints
            for p in range(num_pos):
                for t in range(num_time):
                    save_name = '{0}/{1}_{2}_t{3:05d}xy{4:03d}c3.tif'.format(
                        d, DATE, bn, t, p)

                # Save the mean image.
                skimage.io.imsave(save_name, mean_im)
