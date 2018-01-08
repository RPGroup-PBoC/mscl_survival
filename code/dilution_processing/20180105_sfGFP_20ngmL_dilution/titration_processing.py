import numpy as np
import skimage.io
import skimage.morphology
import skimage.measure
import matplotlib.pyplot as plt
import mscl.plotting
import mscl.image
import pandas as pd
import glob
colors = mscl.plotting.set_plotting_style()

# Load the images.
DATE = 20171207
BASE_STRAIN = 'sfGFP_20ngmL'
IP_DIST = 0.16
data_dir = '../../../data/images/{0}_{1}_dilution/'.format(
    DATE, BASE_STRAIN)

# Get all of the BF and GFP files.
bf_files = np.sort(
    glob.glob('{0}/titration/*/*/*Brightfield*.tif'.format(data_dir)))
gfp_files = np.sort(glob.glob('{0}/titration/*/*/*GFP*.tif'.format(data_dir)))

# Set up the DatFrame.
columns = ['date', 'strain', 'atc_conc', 'mean_int', 'area', 'eccentricity']
df = pd.DataFrame([], columns=columns)

for bf, gfp in zip(bf_files, gfp_files):
    # Get the identifying information.
    try:
        date, atc_conc, _, _ = bf.split('/')[-3].split('_')
        strain = 'dilution'
        atc_conc = int(atc_conc.split('ngmL')[0])
    except:
        date, strain, _ = bf.split('/')[-3].split('_')
        atc_conc = 0

        # Load the two images.
    bf_im = skimage.io.imread(bf)
    gfp_im = skimage.io.imread(gfp)

    # Segment the phase image using otsu.
    seg = mscl.image.threshold_phase(bf_im)

    # Extract the properties.
    props = skimage.measure.regionprops(seg, gfp_im)
    for p in props:
        area = p.area * IP_DIST**2
        ecc = p.eccentricity
        mean_int = p.mean_intensity

        # Assemble the dict and add it to the dataframe.
        cell_dict = dict(date=date, strain=strain, atc_conc=atc_conc,
                         mean_int=mean_int, area=area, eccentricity=ecc)
        df = df.append(cell_dict, ignore_index=True)

# %%
# Add the delta measurements.
delta_files = glob.glob('{0}/photobleaching/*delta*'.format(data_dir))
delta_files

# %%
# Look at the distributions of intensity.
mean_auto = df.loc[df['strain'] == 'autofluorescence']['mean_int'].mean()
dilution_df = df.loc[df['strain'] == 'dilution']
dilution_df['sub'] = dilution_df['mean_int'] - mean_auto
grouped = dilution_df.groupby(['atc_conc'])[['atc_conc', 'sub']].mean()
grouped

# Plot the change in mean intensity with atc concentration for the dilution strain.
fig, ax = plt.subplots(1, 1)
ax.plot(grouped)
