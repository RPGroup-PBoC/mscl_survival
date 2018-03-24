# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import skimage.io
import skimage.morphology
import skimage.segmentation
import pandas as pd
import glob
import sys
sys.path.insert(0, '../')
import mscl.plotting
import mscl.stats
import mscl.image
import mscl.process
colors = mscl.plotting.set_plotting_style()
FIG_NO = 6
IP_DIST = 0.16
%matplotlib inline
files = glob.glob('../data/images/*/Pos*/*GFP*.tif')
fluo_ims = skimage.io.ImageCollection(files, conserve_memory=False)

# Load speciic images and flatten.
bf_im = skimage.io.imread(
    glob.glob('../data/images/*/Pos3/*Brightfield*.tif')[0])
gfp_im = skimage.io.imread(glob.glob('../data/images/*/Pos3/*GFP*.tif')[0])
flat_im = mscl.image.generate_flatfield(gfp_im, np.median(fluo_ims, axis=0))

fig, ax = plt.subplots(1, 1)
ax.imshow(flat_im, cmap='inferno')
ax.set_axis_off()
ax.hlines(-6, 0, 10 / IP_DIST, 'k')
plt.savefig('figS{}_flattened_gfp_image.svg'.format(
    FIG_NO), bbox_inches='tight')

fig, ax = plt.subplots(1, 1)
ax.imshow(gfp_im, cmap='inferno')
ax.set_axis_off()
ax.hlines(-6, 0, 10 / IP_DIST, 'k')
plt.savefig('figS{}_raw_gfp_image.svg'.format(FIG_NO), bbox_inches='tight')


fig, ax = plt.subplots(1, 1)
ax.imshow(np.median(fluo_ims, axis=0), cmap='inferno')
ax.set_axis_off()
ax.hlines(-6, 0, 10 / IP_DIST, 'k')
plt.savefig('figS{}_illumination_profile_image.svg'.format(
    FIG_NO), bbox_inches='tight')
# ---- Segmentation -------------------
# Generate the segmentation
seg = mscl.image.contour_seg(bf_im, level=0.3, min_int=0.4)
seg_boundaries = skimage.measure.find_contours(seg.astype('bool'), level=0)

# Save the example segmentation.
fig, ax = plt.subplots(1, 1)
ax.imshow(bf_im, cmap='Greys_r')

color_list = sns.color_palette('husl', n_colors=np.max(seg) + 1)
for i, c in enumerate(seg_boundaries):
    ax.plot(c[:, 1], c[:, 0], color=color_list[i], lw=1)
    ax.fill(c[:, 1], c[:, 0], color=color_list[i], alpha=0.5)
ax.hlines(-6, 0, 10 / IP_DIST, 'k')
ax.set_axis_off()
plt.tight_layout()
plt.savefig('figS{}_example_segmentation.svg'.format(
    FIG_NO), bbox_inches='tight')


#%%---- Show Markers -------------------
markers = glob.glob('../data/images/*/Pos3.xml')[0]
ids = mscl.image.marker_parse(markers)

fig, ax = plt.subplots(1, 1)
ax.set_axis_off()

_ = ax.imshow(bf_im, cmap='Greys_r')
for i in range(len(ids)):
    slc = ids.iloc[i]
    if slc['survival'] == True:
        c = 'dodgerblue'
    else:
        c = 'tomato'
    _ = ax.plot(slc['x_pos'], slc['y_pos'], marker='.', ms=8, color=c)

_ = ax.hlines(-6, 0, 10 / IP_DIST, 'k')
plt.tight_layout()
plt.savefig('figS{}_markers.svg'.format(FIG_NO), bbox_inches='tight')


#%% --- Map to Fluorescence --------------
link_df = mscl.image.link_markers(ids, seg, gfp_im)

# Mask the fluorescence image.
blank_im = np.zeros_like(flat_im)
blank_im += (seg.astype(bool) * flat_im)

# Instantiate the figure and format
fig, ax = plt.subplots(1, 1)
ax.set_axis_off()
_ = ax.hlines(-6, 0, 10 / IP_DIST, 'k')
_ = ax.imshow(blank_im, vmin=np.min(flat_im[seg.astype(bool)]),
              vmax=0.95 * np.max(flat_im[seg.astype(bool)]),
              cmap='inferno')

# Selectively color the segmentation.
color_dict = {True: 'dodgerblue', False: 'tomato'}
for i, b in enumerate(seg_boundaries):
    cell = link_df[link_df['mask_label'] == i + 1]
    id = cell['survival'].values[0]
    c = color_dict[id]
    _ = ax.plot(b[:, 1], b[:, 0], lw=1, color=c)

plt.tight_layout()
plt.savefig('figS{}_linked_markers.svg'.format(FIG_NO), bbox_inches='tight')


mscl.image.link_markers(ids, )
link_df['nc'] = 2 * link_df['intensity'] * 5.78 / 4400
link_df
