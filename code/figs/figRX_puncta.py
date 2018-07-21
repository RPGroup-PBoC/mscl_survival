# %%
# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
import skimage.io
import skimage.measure
import skimage.morphology
sys.path.insert(0, '../../')
import mscl.image
import mscl.stats
import mscl.plotting
colors = mscl.plotting.set_plotting_style()


# Process the high expression data. 
files = glob.glob('../../data/images/20170421_sd1_pre_*/*.ome.tif')
ims = skimage.io.ImageCollection(files)

# Generate a median image to flatten.
fluos = [ims[i][1] for i in range(len(ims))]
phase = [ims[i][0] for i in range(len(ims))]
med_fluo = np.median(fluos, axis=0)

# Flatten each fluorescent image.
fluo_flat = [mscl.image.generate_flatfield(im, med_fluo) * 4 for im in fluos]

# %%
# Loop through each phase image, segment, and collect fluorescence statistics
df = pd.DataFrame([])
sd1_segs = []
for i, im in enumerate(phase):
    seg = mscl.image.contour_seg(phase[0], level=0.2, min_int=0.4)

    # Compute the mean background
    inv_seg = seg < 1
    mean_bg = np.mean(fluo_flat[i][(inv_seg * fluo_flat[i]) > 0])

    # Extract the properties.
    props = skimage.measure.regionprops(seg, fluo_flat[i])
    for p in props:
        var_int = np.var(fluo_flat[i][((seg==p.label) * fluo_flat[i]) > 0] - mean_bg)
        cell_dict = dict(var_int=var_int, mean_int=p.mean_intensity, area=p.area, rbs='sd1',
        cell_label=p.label, file=files[i]) 
        df = df.append(cell_dict, ignore_index=True)
        sd1_segs.append((fluo_flat[i] * (seg==p.label))[p.bbox[0]-20:p.bbox[2]+20, p.bbox[1]-20:p.bbox[3]+20])

# %% Process the low expression samples
files = glob.glob('../../data/images/*sd4_pre_*/*.ome.tif')
ims = skimage.io.ImageCollection(files)

# Generate a median image to flatten.
fluos = [ims[i][1] for i in range(len(ims))]
phase = [ims[i][0] for i in range(len(ims))]
med_fluo = np.median(fluos, axis=0)

# Flatten each fluorescent image.
fluo_flat = [mscl.image.generate_flatfield(im, med_fluo) for im in fluos]
sd4_segs = []
for i, im in enumerate(phase):
    seg = mscl.image.contour_seg(phase[0], level=0.2, min_int=0.4)

    # Compute the mean background
    inv_seg = seg < 1
    mean_bg = np.mean(fluo_flat[i][(inv_seg * fluo_flat[i]) > 0])

    # Extract the properties.
    props = skimage.measure.regionprops(seg, fluo_flat[i])
    for p in props:
        var_int = np.var(fluo_flat[i][((seg==p.label) * fluo_flat[i]) > 0] - mean_bg)
        cell_dict = dict(var_int=var_int, mean_int=p.mean_intensity, area=p.area, rbs='sd4',
        cell_label=p.label, file=files[i]) 
        df = df.append(cell_dict, ignore_index=True)
        sd4_segs.append((fluo_flat[i] * (seg==p.label))[p.bbox[0]-20:p.bbox[2]+20, p.bbox[1]-20:p.bbox[3]+20])

df.to_csv('../../data/csv/sd1_sd4_intensity_variance.csv', index=False)

# %%
# Normalize the mean and variance
normed_dfs = []
grouped = df.groupby('rbs')
for g, d in grouped:
    d = d.copy()
    d['norm_mean'] = d['mean_int'] - d['mean_int'].mean()
    d['noise'] = d['var_int'] / d['mean_int']
    normed_dfs.append(d)
normed_df = pd.concat(normed_dfs)

# %%
fig, ax = plt.subplots(2, 2, figsize=(6, 4))
grouped = normed_df.groupby('rbs')
color_dict = {'sd1':colors['red'], 'sd4':colors['blue']}
axes = {'sd1':ax[0, 0], 'sd4':ax[1, 0]}

bins = np.linspace(-2, 4, 50)
for g, d in grouped:
    axes[g].hist(d['norm_var'], bins=bins, color=color_dict[g], alpha=0.75)

ax[1, 0].set_xlabel('log variance in cell intensity', fontsize=8)
ax[0, 0].set_title(r'SD1 | $\langle N \rangle \approx 500$ [channels / cell]',
                backgroundcolor=colors['pale_yellow'], fontsize=8, y=1.05)
ax[1, 0].set_title(r'SD4 | $\langle N \rangle \approx 80$ [channels / cell]',
                backgroundcolor=colors['pale_yellow'], fontsize=8, y=1.05)
for a in ax.ravel()[::2]:
    a.set_ylabel('counts', fontsize=8)
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
ax[0, 1].axis('off')
ax[1, 1].axis('off')
plt.tight_layout()
plt.savefig('../../figs/figRX_cell_intensity_variance.pdf')



#%% Look at some example cells.
high_ims = []
low_ims = []
high_labels = []
low_labels = []

grouped = normed_df.groupby('rbs')
for g, d in grouped:
    print(g)
    d = d.copy()
    d.sort_values('norm_var', inplace=True)
    high = d.iloc[-6]
    low = d.iloc[5]
    high_ims.append(high['file'])
    high_labels.append(high['cell_label'])
    low_ims.append(low['file'])
    low_labels.append(low['cell_label'])


high_seg = []
low_seg = []
for i, im in enumerate(high_ims):
    im = skimage.io.imread(im)
    phase = im[0, :, :]
    fluo = im[1, :, :]
    seg = mscl.image.contour_seg(phase, level=0.2, min_int=0.4)
    props = skimage.measure.regionprops(seg)
    bbox = props[int(high_labels[i] - 1)].bbox
    high_seg.append(((seg==p.label) * fluo)[bbox[0]-20:bbox[2]+20,  bbox[1]-20:bbox[3]+20])

for i, im in enumerate(low_ims):
    im = skimage.io.imread(im)
    phase = im[0, :, :]
    fluo = im[1, :, :]
    seg = mscl.image.contour_seg(phase, level=0.2, min_int=0.4)
    props = skimage.measure.regionprops(seg)
    bbox = props[int(low_labels[i] - 1)].bbox
    low_seg.append(((seg==p.label) * fluo)[bbox[0]-20:bbox[2]+20,  bbox[1]-20:bbox[3]+20])


fig, ax = plt.subplots(2, 2)
for i, im in enumerate(high_seg):
    if i == 0:
        mult = 4
    else:
        mult =1 
    ax[0, i].imshow(high_seg[i] * mult, vmin=np.min(high_seg[i][high_seg[i] > 0]))
                                     
for i, im in enumerate(low_seg):
    if i == 0:
        mult = 4
    else:
        mult =1 
    ax[1, i].imshow(low_seg[i] * mult, vmin=np.min(low_seg[i][low_seg[i] > 0]))
                                     
# %%
fig, ax = plt.subplots(10, 10, figsize=(6, 6), sharex=True, sharey=True)
ax = ax.ravel()
for i, im in enumerate(sd1_segs):
    ax[i].imshow(im, cmap=plt.cm.viridis)

for a in ax:
    a.grid(False)
plt.tight_layout()
