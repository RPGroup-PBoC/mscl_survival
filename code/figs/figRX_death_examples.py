# -*- coding: utf-8 -*-
# %%
import sys
import numpy as np
import matplotlib.pyplot as plt
import skimage.io
import skimage.morphology
import skimage.segmentation
import skimage.measure
import glob
sys.path.insert(0, '../../')
import mscl.image
import mscl.plotting
colors = mscl.plotting.set_plotting_style()

# Load the image.
files = np.sort(glob.glob('../../data/images/20170518_12SD2_8000ulmin_post/*.tif'))
im = skimage.io.ImageCollection(files)
init_im = im[0]

# Segment the image.
thresh = skimage.filters.threshold_otsu(init_im)
seg  = mscl.image.contour_seg(init_im, level=0.2, min_int=0.4)
props = skimage.measure.regionprops(lab)
# plt.imshow(lab==5)
# %%
# ex_labels = [8, 17, 18, 10, 11,20]
# ex_labels = [3, 4, 5, 6, 7, 8]
ex_labels = [prop.label for prop in props]
# ex_labels = [2, 3, 1, 5, 6, 7, 8, 4]
locs = {loc:i for i, loc in enumerate(ex_labels)}

bbox = []
orientation = []
buffer = 35 
for l in ex_labels: 
    p = props[l-1]
    bb = p.bbox
    bbox.append([int(p.centroid[0] - buffer), int(p.centroid[1] - buffer), 
                    int(p.centroid[0] + buffer), int(p.centroid[1] + buffer)])
    orientation.append(p.orientation)


# Instantiate figure and format axes
fig, ax = plt.subplots(len(ex_labels), 10, sharex=True, 
                    sharey=True, figsize=(6, 5))
ax[0, 4].set_title(r'12SD2 $\langle N \rangle \approx 300$ effective channels per cell  | fast shock (> 1.0 Hz)', fontsize=8, backgroundcolor=colors['pale_yellow'],
                    y=1.18)
for a in ax.ravel():
    a.set_xticks([])                
    a.set_yticks([])                

# Add appropriate labels
fig.text(0.45, 0.15, 'time [min]', fontsize=8)

for i in range(len(ex_labels)):
    bb = bbox[i]
    c = 0
    for j in range(0, len(files), round(len(files) / 10)):
        pos = im[j]
        ax[locs[ex_labels[i]], c].imshow(pos[bb[0]:bb[2], bb[1]:bb[3]], cmap=plt.cm.Greys_r,
                        vmin=3000, vmax=12000)
        if i == len(ex_labels) - 1:
            ax[-1, c].set_xlabel(j, fontsize=8) 
        c += 1
plt.subplots_adjust(hspace=-.62, wspace=0.06)
plt.savefig('../../figs/figRX_death_example.pdf')




