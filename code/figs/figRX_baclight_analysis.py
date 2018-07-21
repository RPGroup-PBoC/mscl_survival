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
segs = []
fl = []
ims = skimage.io.ImageCollection(pos)
baclight = skimage.io.imread(bl_pos[0])
seg = mscl.image.contour_seg(ims[0], level=0.2, min_int=0.4)
props = skimage.measure.regionprops(seg)
for p in props:
    x, y = p.centroid
    bb = np.s_[int(x-50):int(x+50), int(y-50):int(y+50)]
    for i in ims:
        segs.append(i[bb])
    fl.append(baclight[bb]) 

