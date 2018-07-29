# -*- coding: utf-8 -*-
# %%
import sys
import numpy as np
import matplotlib.pyplot as plt
import skimage.io
import skimage.morphology
import skimage.measure
import glob
import pandas as pd
import tqdm
sys.path.insert(0, '../../')
import mscl.stats
import mscl.image
import mscl.plotting
colors = mscl.plotting.set_plotting_style()

# Define the data directory.
data_dir = '../../data/images/20180719_baclight_staining/'

#Grab all positions.
pos = np.sort(glob.glob('{}Pos12/*.tif'.format(data_dir)))
bl_pos = np.sort(glob.glob('{}baclight/Pos12/*Bac*.tif'.format(data_dir)))
ids = {}
fl = []
ims = skimage.io.ImageCollection(pos)
baclight = skimage.io.imread(bl_pos[0])
seg = mscl.image.contour_seg(ims[0], level=0.2, min_int=0.4)
props = skimage.measure.regionprops(seg)
labels = [p.label for p in props]
fig, ax = plt.subplots(8, 12, figsize=(6, 6), sharex=True, sharey=True)
for a in ax.ravel():
    a.grid(False)
    a.set_xticks([])
    a.set_yticks([])

for i in range(len(labels)):
    x, y = props[labels[i] - 1].centroid
    bb = np.s_[int(x-30):int(x+30), int(y-30):int(y+30)]
    for j, im in enumerate(ims[::12]):
        ax[i-1, j].imshow(im[bb], cmap=plt.cm.Greys_r, vmin=np.min(ims[0]), vmax=np.max(ims[0]))
        ax[-1, j].set_xlabel(j * 12, fontsize=8)
    ax[i-1, -1].imshow(baclight[bb], cmap='viridis', vmin=np.min(baclight), vmax=np.max(baclight))   
    
plt.subplots_adjust(hspace=-.79, wspace=0.1)
fig.text(0.45, 0.2, 'time [min]', fontsize=8)
plt.savefig('../../figs/figRX_dye_example.pdf', bbox_inches='tight')