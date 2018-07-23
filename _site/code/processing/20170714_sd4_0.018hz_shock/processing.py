import numpy as np
import glob
import skimage.io
import mscl_utils as mscl
import bokeh.io
import bokeh.plotting
import pandas as pd
import os
import tqdm


# %% Define the relevant parameters
DATE = 20170714
RBS = 'sd4'
FLOW_RATE = '0.018'  # in units of Hz.
EXP = 'shock'
run_number = 1
# EXP = 'cal'
# REP = [0000, 0001]

# --------------------------------------------------------------------------
# Shouldn't need to touch information below this point
# --------------------------------------------------------------------------
ip_dist = 0.160  # in units of µm per pixel.
# Automatically parse the directory.
if EXP.lower() == 'shock':
    data_dir = '../../../data/images/{0}_{1}_pre_*_{2}hz_shock/'.format(
        DATE, RBS, FLOW_RATE)
elif EXP.lower() == 'cal':
    if type(REP) is not list:
        data_dir = '../../../data/images/{0}_{1}_pre_{2}_intensity_calibration'.format(
            DATE, RBS, REP)
    else:
        data_dir = '../../../data/images/{0}_{1}_pre_*_intesity_calibration'.format(
            DATE, RBS)

# Grab the files and process.
data_dir
tif_files = glob.glob(data_dir + '/*.ome.tif')
if len(tif_files) > 0:
    _ims = skimage.io.ImageCollection(tif_files, conserve_memory=False)
    phase_ims = [i[0] for i in _ims]
    fluo_ims = [i[1] for i in _ims]
    metadata = glob.glob(data_dir + '/*.txt')
else:
    phase_files = glob.glob(data_dir + '/Pos*/*Brightfield*.tif')
    phase_ims = skimage.io.ImageCollection(phase_files, conserve_memory=False)
    gfp_files = glob.glob(data_dir + '/Pos*/*GFP*.tif')
    fluo_ims = skimage.io.ImageCollection(gfp_files, conserve_memory=False)
    metadata = glob.glob(data_dir + '/Pos*/*.txt')

# Load the mareker files if it's a shock experiment.
if EXP.lower() == 'shock':
    markers = glob.glob(data_dir + '/*.xml')

# Generate the averiage illumination profile.
ff_ims = mscl.median_flatfield(fluo_ims, medfilter=True)

# Choose the random segmentation file.
rand = np.random.choice(np.arange(0, len(phase_ims), 1))
dfs = []
for i, ph in tqdm.tqdm(enumerate(phase_ims), desc='processing images',
                       total=len(phase_ims)):
    seg = mscl.contour_seg(ph, level=0.2, min_int=0.3)
    seg, num_obj = skimage.measure.label(seg > 0, return_num=True)

    # Determine if the segmentation should be saved.
    if i == rand:
        mscl.save_seg(
            'output/{0}_{1}_example_segmentation.png'.format(DATE, RBS), ph, seg)

    if num_obj > 0:
        # Compute the background information
        mean_bg = mscl.compute_mean_bg(ph, ff_ims[i])
        mean_bg = mean_bg / ip_dist**2
        exposure = mscl.scrape_metadata(metadata[i], return_date=False)

        # For shock experiments, link the markers.
        if EXP.lower() == 'shock':
            _markers = mscl.marker_parse(markers[i])
            _df = mscl.link_markers(_markers, seg, ff_ims[i])
            if i == rand:
                mscl.show_connections(
                    'output/{0}_{1}_example_connections.png'.format(DATE, RBS), ph, _df)
            _df['flow_rate'] = float(FLOW_RATE)
            _df.drop(['x_pos', 'y_pos', 'mask_label', 'label_cent_x',
                      'label_cent_y', 'dist'], axis=1, inplace=True)

        else:
            props = skimage.measure.regionprops(seg, ff_Ims[i])
            intensity = [prop.mean_intensity / ip_dist**2 for prop in props]
            area = [prop.area * ip_dist**2 for prop in props]
            _df = pd.DataFrame(dict(intensity=intensity, area=area))

        # Insert the other necessary information.
        _df['date'] = DATE
        _df['rbs'] = RBS
        _df['mean_bg'] = mean_bg
        _df['exposure_ms'] = exposure['GFP_exp_ms']

        # Append the dataframe to the list
        dfs.append(_df)

df = pd.concat(dfs, axis=0)

# Save the dataframe to disk.
if EXP.lower() == 'shock':
    target = 'output/{0}_{1}_r{2}_{3}Hz_shock.csv'.format(DATE, RBS,
                                                          run_number,
                                                          FLOW_RATE)
else:
    target = 'output/{0}_{1}_r{2}_intensity_calibration.csv'.format(DATE, RBS,
                                                                    run_number)

if os.path.exists(target) is True:
    os.remove(target)
with open('comments.txt', 'r') as f:
    comments = f.readlines()
with open(target, 'a') as f:
    for line in comments:
        f.write(line)
    df.to_csv(f, index=False)


#%% Plot the Histograms and ECDFs of survival.
p1 = mscl.bokeh_boiler(
    **dict(y_axis_label='ECDF', x_axis_label='raw intensity (a.u.)',
           x_axis_type='log', y_range=[0, 1]))
p2 = mscl.bokeh_boiler(
    **dict(y_axis_label='count', x_axis_label='raw intensity (a.u.)',
           x_axis_type='log'))
n_bins = np.max([int(np.sqrt(len(df))), 25])
bins = np.linspace(df['intensity'].min(), df['intensity'].max(), n_bins)
if EXP.lower() == 'shock':
    grouped = df.groupby('survival')
    colors = {'True': 'slategray', 'False': 'tomato'}
    labels = {'True': 'survival', 'False': 'death'}
    for g, d in grouped:
        x = np.sort(d['intensity'])
        y = np.arange(0, len(d)) / len(d)
        hist, edges = np.histogram(d['intensity'], bins=bins)
        p1.line(x, y, color=colors[str(g)],
                legend=labels[str(g)], line_width=2)
        p1.circle(x, y, color=colors[str(g)], legend=labels[str(g)],
                  size=5)

        p2.quad(bottom=0, top=hist,
                left=edges[:-1], right=edges[1:], color=colors[str(g)], alpha=0.5)

    p1.legend.background_fill_color = '#FFEDCE'
    p1.legend.location = 'bottom_right'
else:
    hist, edges = np.histogram(df['intensity'], bins=bins)
    x = np.sort(df['intensity'])
    y = np.arange(0, len(df['intensity']), 1) / len(df['intensity'])
    p1.line(x, y, color='dodgerblue', line_width=2)
    p1.circle(x, y, color='dodgerblue')
    p2.quad(bottom=0, top=hist, left=edges[:-1], right=edges[1:])

p1.title.text = '{0} {1} intensity distributions'.format(DATE, RBS)
layout = bokeh.layouts.row(p1, p2)
fname = 'output/{0}_{1}_intensity_distribution.'.format(DATE, RBS)
bokeh.io.output_file(fname + 'html')
bokeh.io.save(layout)
p1.toolbar_location = None
p2.toolbar_location = None
bokeh.io.export_png(layout, filename=fname + 'png')
