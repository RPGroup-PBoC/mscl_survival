# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
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
gfp_im = skimage.io.imread(glob.glob('../data/images/*/Pos0/*GFP*.tif')[0])
flat_im = mscl.image.generate_flatfield(gfp_im, np.median(fluo_ims, axis=0))


# ---- Segmentation -------------------
# Generate the segmentation
seg = mscl.image.contour_seg(bf_im, level=0.2, min_int=0.3)
seg_boundaries = skimage.measure.find_contours(seg.astype('bool'), level=0)

# Save the example segmentation.
fig, ax = plt.subplots(1, 1)
ax.imshow(bf_im, cmap='Greys_r')
color_list = (list(colors.values()) + list(colors.values()) +
              list(colors.values()))[::3]
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
new_seg = np.zeros_like(seg)
for _, id in enumerate(link_df['mask_label'].values):
    new_seg += (seg == id)

plt.imshow(new_seg)

blank_im = np.zeros_like(gfp_im)
blank_im += (seg.astype(bool) * gfp_im)

# Instantiate the figure and format
fig, ax = plt.subplots(1, 1)
ax.set_axis_off()
_ = ax.vlines(-6, 0, 10 / IP_DIST, 'k')
_ = ax.imshow(blank_im, vmin=np.min(gfp_im[seg.astype(bool)]),
              cmap='inferno')

# Selectively color the segmentation.
color_dict = {True: 'dodgerblue', False: 'tomato'}
for i, b in enumerate(seg_boundaries):
    id = link_df[link_df['mask_label'] == i + 1]['survival'].values[0]
    c = color_dict[id]
    _ = ax.plot(b[:, 1], b[:, 0], lw=1, color=c)
