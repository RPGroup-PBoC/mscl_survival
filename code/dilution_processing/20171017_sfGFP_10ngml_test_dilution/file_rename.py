import glob
import os
import shutil

# Define the experiment parameters.
DATE = 20171017
BASENAME = 'sfGFP_10ngmL'
channels = ['Brightfield', 'GFP']
channel_dict = {'Brightfield': 1, 'GFP': 2}
keep_origs = True

# ----------------------------------------------------------------------------
# Nothing below here should change.
# ----------------------------------------------------------------------------
data_dir = '../../../data/images/{0}_{1}_dilution/'.format(DATE, BASENAME)

# Scrape positions.
positions = glob.glob('{0}Pos*/'.format(data_dir))

# Loop through positions and grab file names.
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
positions = glob.glob('{0}{1}_autofluorescence*/Pos*'.format(data_dir, DATE))

if os.path.isdir('{0}autofluorescence'.format(data_dir)) is False:
    os.mkdir('{0}autofluorescence'.format(data_dir))

for p in positions:
    # Parse the position.
    xy = '{0:03d}'.format(int(p.split('/')[-1].split('s')[1]))
    for ch in channels:
        files = glob.glob(p + '/*{0}*.tif'.format(ch))
        files.sort()
        for f in files:
            # Parse the relevant info
            _, t, c, _ = f.split('/')[-1].split('_')

            # Define the new name and copy.
            new_name = '{0}_{1}_autofluorescence_t{2:05d}xy{3:03d}c{4}.tif'.format(
                DATE, BASENAME, int(t), int(xy), channel_dict[ch])

            shutil.copy(
                f, '{0}autofluorescence/{1}'.format(data_dir, new_name))

if os.path.isdir('{0}autofluorescence/originals'.format(data_dir)) is False:
    os.mkdir('{0}autofluorescence/originals'.format(data_dir))

for p in positions:
    shutil.move(p, '{0}autofluorescence/originals/{1}'.format(data_dir,
                                                              p.split('/')[-1]))

# Determine if originals should be kept.
if keep_origs is False:
    shutil.rmtree('{0}autoflurescence/originals'.format(data_dir))
