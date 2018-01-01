import glob
import os
import shutil
import skimage.io
import mscl_utils as mscl
import skimage.io
import numpy as np

# Define the experiment parameters.
DATE = 20171206
BASENAME = 'sfGFP_10ngmL'
channels = ['Brightfield', 'GFP']
channel_dict = {'Brightfield': 1, 'GFP': 2}
keep_origs = True
generate_flatfield = True
# ----------------------------------------------------------------------------
# Nothing below here should change.
# ----------------------------------------------------------------------------
data_dir = '../../../data/images/{0}_{1}_dilution/'.format(DATE, BASENAME)

# Assemble the mean illumination image if appropriate.
if generate_flatfield is True:
    # Load the slide images.
    slide_ims = glob.glob(
        '{0}{1}_fluorescent_slide*/Pos*/*.tif'.format(data_dir, DATE))
    ims = skimage.io.ImageCollection(slide_ims)
    mean_im = mscl.average_stack(ims)
    flatfield = {'Brightfield': False,  'GFP': True}
else:
    flatfield = {'Brightfield': False, 'GFP': False}

# Scrape positions.
positions = glob.glob('{0}{1}_growth*/Pos*/'.format(data_dir, DATE))
for p in positions:
    # Parse the position.
    xy = '{0:03d}'.format(int(p.split('/')[-2].split('s')[1]))

    # Loop through provided channels and rename images.
    for ch in channels:
        files = glob.glob(p + '*{0}*.tif'.format(ch))
        files.sort()

        for f in files:
            # Parse the relevant info
            _, t, c, _ = f.split('/')[-1].split('_')

            # Define the new name and copy.
            new_name = '{0}_{1}_t{2:05d}xy{3:03d}c{4}.tif'.format(
                DATE, BASENAME, int(t), int(xy), channel_dict[ch])
            if flatfield[ch] is True:
                im = skimage.io.imread(f)
                ff_im = mscl.generate_flatfield(im, mean_im,
                                                median_filt=False)
                skimage.io.imsave('{0}{1}'.format(data_dir, new_name), ff_im)
            else:
                shutil.copy(f, '{0}{1}'.format(data_dir, new_name))

    # Rename and save the metadata files.
    mfiles = glob.glob('{0}metadata.txt'.format(p))
    mfiles
    for m in mfiles:
        new_m_name = 'xy{0:03d}_metadata.txt'.format(int(xy))
        shutil.copy(m, '{0}{1}'.format(data_dir, new_m_name))

# Move all original files into a separate folder.
if os.path.isdir('{0}/originals'.format(data_dir)) is False:
    os.mkdir('{0}/originals'.format(data_dir))

for p in positions:
    shutil.move(p, '{0}{1}/{2}'.format(data_dir,
                                       'originals', p.split('/')[-2]))

# Determine if originals should be kept.
if keep_origs is False:
    shutil.rmtree('{0}{1}'.format(data_dir, 'originals'))


# %% Do the same for the autofluroescence files.
positions = glob.glob('{0}{1}_autofluorescence*/*.tif'.format(data_dir, DATE))

if os.path.isdir('{0}autofluorescence'.format(data_dir)) is False:
    os.mkdir('{0}autofluorescence'.format(data_dir))

for p in positions:
    # Parse the position.
    xy = '{0:02d}'.format(int(p.split('/')[-1].split('s')[-1].split('.')[0]))

    # Define the new names.
    new_bf_name = '{0}_{1}_autofluorescence_t{2:05d}xy{3:03d}c{4}.tif'.format(DATE, BASENAME,
                                                                              0, int(
                                                                                  xy),
                                                                              'Brightfield')
    new_fl_name = '{0}_{1}_autofluorescence_t{2:05d}xy{3:03d}c{4}.tif'.format(DATE, BASENAME,
                                                                              0, int(
                                                                                  xy),
                                                                              'GFP')
    # Load the image.
    im = skimage.io.imread(p)

    # Split into phase and gfp channels.
    bf_im = im[0, :, :]
    fl_im = im[1, :, :]
    np.shape(im)

    if flatfield['GFP'] is True:
        ff_im = mscl.generate_flatfield(fl_im, mean_im, median_filt=False)
        skimage.io.imsave('{0}autofluorescence/{1}'.format(data_dir, new_fl_name,),
                          ff_im)
    else:
        shutil.copy(
            fl_im, '{0}autofluorescence/{1}'.format(data_dir, new_fl_name))

    # Rename the brightfield image.
    shutil.copy(bf_im, '{0}autofluorescence/{1}'.format(data_dir, new_bf_name))


if os.path.isdir('{0}autofluorescence/originals'.format(data_dir)) is False:
    os.mkdir('{0}autofluorescence/originals'.format(data_dir))

for p in positions:
    shutil.move(p, '{0}autofluorescence/originals/{1}'.format(data_dir,
                                                              p.split('/')[-1]))

# Determine if originals should be kept.
if keep_origs is False:
    shutil.rmtree('{0}autoflurescence/originals'.format(data_dir))
