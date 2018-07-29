# %%
# -*- coding: utf-8 -*-
import sys
import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
import scipy.stats
import skimage.io
import skimage.measure
import skimage.morphology
import scipy.stats
sys.path.insert(0, '../../')
import mscl.image
import mscl.stats
import mscl.plotting
colors = mscl.plotting.set_plotting_style()



# %%
# Work out the Binomial distribution for a patch of 250 x 250 nm 
sheet_size = 4 # in µm sq
box_size = 0.250**2 # in µm sq
n_sites = sheet_size / box_size
p = 1 / n_sites
N = 500
chan_range = np.arange(0, N)
pmf = scipy.stats.binom.pmf(k=chan_range, n=N, p=p)
np.sum(pmf[12:])
# Istantiate the figure and format axes
fig, ax = plt.subplots(1, 1, figsize=(2.3, 2.1))
ax.set_xlim([0, 30])
ax.set_xlabel('MscL channels per site', fontsize=8)
ax.set_ylabel('probability', fontsize=8)
ax.xaxis.set_tick_params(labelsize=8)
ax.yaxis.set_tick_params(labelsize=8)

# plot the pmf
ax.plot(chan_range, pmf, '.', color=colors['red'], ms=4)
for c, p in zip(chan_range, pmf):
    ax.vlines(c, 0, p, lw=1, color=colors['red'])
plt.tight_layout()
plt.savefig('../figs/figRX_pmf_channels.pdf')




# %%
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
# Compute the intensity of puncta.

# %%
# Loop through each phase image, segment, and collect fluorescence statistics
df = pd.DataFrame([])
sd1_segs = []
high_subims = []
for i, im in enumerate(phase):
    seg = mscl.image.contour_seg(phase[0], level=0.2, min_int=0.4)

    punct_seg = skimage.morphology

    # Compute the mean background
    inv_seg = seg < 1
    mean_bg = np.mean(fluo_flat[i][(inv_seg * fluo_flat[i]) > 0])

    # Extract the properties.
    props = skimage.measure.regionprops(seg, fluo_flat[i])
    for p in props:
        bb = np.s_[p.bbox[0]-20:p.bbox[2]+20, p.bbox[1]-20:p.bbox[3]+20]
        high_subims.append(((seg==p.label) * (fluo_flat[i] - mean_bg))[bb] * 4) 
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
low_subims = []
sd4_segs = []
for i, im in enumerate(phase):
    seg = mscl.image.contour_seg(phase[0], level=0.2, min_int=0.4)

    # Compute the mean background
    inv_seg = seg < 1
    mean_bg = np.mean(fluo_flat[i][(inv_seg * fluo_flat[i]) > 0])
    
    # Extract the properties.
    props = skimage.measure.regionprops(seg, fluo_flat[i])
    for p in props:
        bb = np.s_[p.bbox[0]-20:p.bbox[2]+20, p.bbox[1]-20:p.bbox[3]+20]
        low_subims.append((fluo_flat[i] - mean_bg)[bb] * 4) 
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
    d['norm_var'] = d['var_int'] / d['mean_int']
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
    high = d.iloc[-1]
    low = d.iloc[2]
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
    high_seg.append(fluo[bbox[0]-20:bbox[2]+20,  bbox[1]-20:bbox[3]+20])

for i, im in enumerate(low_ims):
    im = skimage.io.imread(im)
    phase = im[0, :, :]
    fluo = im[1, :, :]
    seg = mscl.image.contour_seg(phase, level=0.2, min_int=0.4)
    props = skimage.measure.regionprops(seg)
    bbox = props[int(low_labels[i] - 1)].bbox
    low_seg.append(fluo[bbox[0]-20:bbox[2]+20,  bbox[1]-20:bbox[3]+20])


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

# %%
grouped = normed_df.groupby('rbs')
perc = {}
for g, d in grouped:
    high_var = np.sum(d['norm_var'] > 100)
    total = len(d)
    print(g, total)
    perc[g] = high_var/total

# %%
tot = []
for i, _ in enumerate(high_subims):
    # Extract the intensity of segmented cells.
    slc = high_subims[i][high_subims[i] > 0].flatten()
    for s in slc:
        tot.append(s)
fig, ax = plt.subplots(1, 1)
ax.hist(tot, bins=500)
# ax.set_yscale('log')

# %%
# Define the pixel correction factor.
correction_factor = 250**2 / 160**2
ind =0 
channels = []
mean_subims = [] 
perc_per_punctum = []
thresh = 6000 # Determined empirically
n_punct = 0
for i, im in enumerate(high_subims):
    mean_subims.append(np.sum(im[im > 0])/ alpha)
    seg = (im > 6000)
    puncta_seg = scipy.ndimage.binary_fill_holes(seg)
    large = skimage.morphology.remove_small_holes(puncta_seg)
    large = skimage.morphology.remove_small_objects(large, min_size=3)
    structure = skimage.morphology.square(2)
    dil = skimage.morphology.binary_dilation(large, structure)
    label, obj = skimage.measure.label(dil, return_num=True)
    if obj > 0:
        n_punct += obj
    props = skimage.measure.regionprops(label, im * 4)
    for p in props:
        channels.append((p.mean_intensity / alpha) * correction_factor)
        perc_per_punctum.append((p.mean_intensity / alpha) * correction_factor / (np.sum(im[im > 0]) * alpha**-1))

# %%
fig, ax = plt.subplots(1, 2, figsize=(1.8, 2.1))
ax[1].set_ylabel('channel fraction per punctum', fontsize=8)
ax[1].yaxis.tick_right()
ax[1].yaxis.set_label_position('right')
ax[0].set_ylabel('channels per punctum', fontsize=8)
for a in ax:
    a.yaxis.set_tick_params(labelsize=8)

# Add the channels per punctum.
sns.swarmplot(y=channels, ax=ax[0], color=colors['red'])
sns.boxplot(channels, ax=ax[0], orient='v', color=colors['pale_yellow'], linewidth=1)
sns.swarmplot(y=perc_per_punctum, ax=ax[1], color=colors['red'])
sns.boxplot(perc_per_punctum, ax=ax[1], orient='v',  color=colors['pale_yellow'], linewidth=1)
plt.savefig('../figs/figRX_fraction_per_punctum.pdf', bbox_inches='tight')

np.sum(pmf[13:])
n_punct /len(high_subims)
colors.keys()