# -*- coding: utf-8 -*-
# %%
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import  matplotlib.gridspec as gridspec
sys.path.insert(0, '../../')
import mscl.plotting
colors = mscl.plotting.set_plotting_style()

# Load the data.
data = pd.read_csv('../../data/csv/MLG910_12118002.csv')
data.columns = ['time', 'pa', 'mmHg']

# Instantiate the figure.
fig = plt.figure(figsize=(4, 4))
gs = gridspec.GridSpec(3, 1)
ax1 = fig.add_subplot(gs[:2, 0])
ax2 = fig.add_subplot(gs[2, 0])
ax1.set_xlim([19, 23.8])
ax2.set_xlim([19, 23.8])
ax1.set_ylim([0, 525])
ax2.set_ylim([-350, 0])

# Format the axes
ax1.xaxis.set_ticklabels([])
ax1.xaxis.set_tick_params(labelsize=8)
ax2.xaxis.set_tick_params(labelsize=8)
ax1.yaxis.set_tick_params(labelsize=8)
ax2.yaxis.set_tick_params(labelsize=8)

# Add the appropriate labels
ax1.set_ylabel('current [pA]', fontsize=8)
ax2.set_ylabel('pressure\n [mmHg]', fontsize=8)
ax2.set_xlabel('time [s]', fontsize=8)

# Add marker labels
ax1.text(0.08, 0.93, 'MscS', fontsize=8, transform=ax1.transAxes,
        backgroundcolor=colors['pale_blue'])
ax2.text(0.46, 0.93, 'MscL-sfGFP', fontsize=8, transform=ax1.transAxes,
        backgroundcolor=colors['pale_yellow'])

# Plot the traces and color red
_ = ax1.plot(data['time'], data['pa'], '-', color=colors['red'], lw=0.5)
_ = ax2.plot(data['time'], data['mmHg'], '-', color=colors['red'], lw=1)

# Label the MscS points
_ = ax1.vlines(19.6, -1.5, 525, lw=28, color=colors['pale_blue'], zorder=-1)
_ = ax2.vlines(19.6, 4, -350, lw=28, color=colors['pale_blue'], zorder=-1)

# Label the MscL points
_ = ax1.vlines(21.7, -1.5, 525, lw=100, color=colors['pale_yellow'], zorder=-1)
_ = ax2.vlines(21.7, 4, -350, lw=100, color=colors['pale_yellow'], zorder=-1)
plt.savefig('../../figs/figRX_electrophysiology_trace.pdf', bbox_inches='tight')
plt.savefig('../../figs/figRX_electrophysiology_trace.png', bbox_inches='tight', 
            dpi=300)
# %% Compute the conductance
R = 20 # in mV
slc = data[(data['time'] > 20.5) & (data['time'] < 22.8)]
plt.plot(slc['time'], slc['pa'], lw=0.3)
I_open = np.mean(slc[slc['pa'] > 260]['pa'])
I_closed = np.mean(slc[slc['pa'] < 200]['pa'])
conductance = (I_open - I_closed) / R
conductance
I_open - I_closed
