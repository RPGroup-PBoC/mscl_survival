# -*- coding: utf-8 -*-
#%%
import sys
import skimage.io
import skimage.measure
import matplotlib.pyplot as plt
sys.path.insert(0, '../../')
import mscl.plotting
import mscl.image
import glob
colors = mscl.plotting.set_plotting_style()

# Load all mlg910 images.
bf_files = glob.glob('../../data/images/mlg910/*c1*.tif') 
bf_ims = skimage.io.ImageCollection(bf_files)
gfp_files = glob.glob('../../data/images/mlg910/*c2*.tif')
gfp_ims = skimage.io.ImageCollection(gfp_files)
plt.imshow(gfp_ims[0], cmap='magma')

#%% Segment.
seg = mscl.image.contour_seg(bf_ims[0], level=0.3, min_int=1.0)

# Get the properties.
props = skimage.measure.regionprops(seg, gfp_ims[0])


    
plt.imshow(seg)

