# -*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pystan
import sys
sys.path.insert(0, '../')
import mscl.plotting
import mscl.stats
import mscl.mcmc

colors = mscl.plotting.set_plotting_style()

FIG_NO = 'X4'

data = pd.read_csv('../data/csv/compiled_data.csv')
shock_data = data[data['class'] == 'shock'].copy()

#%% -- Bivariate regression --------------
# Compile the stan model.
sm = pystan.StanModel('../code/bivariate_logistic_regression.stan')
data_dict = {'N': len(shock_data), 'predictor_1': shock_data['chan_per_cell'],
             'predictor_2': shock_data['flow_rate'],
             'output': shock_data['survival'].astype(int)}
# Sample and convert the traces to a dataframe
bivariate_traces = sm.sampling(data=data_dict, iter=5000, chains=4)
bivariate_df = mscl.mcmc.chains_to_dataframe(bivariate_traces)

# Compute the means for all coefficients.
_, beta_2, beta_1, beta_0 = np.mean(bivariate_df).values

#%% Set up the grid
chan_range = np.linspace(0, 800, 1000)
shock_range = np.linspace(0, 3, 1000)
X, Y = np.meshgrid(chan_range, shock_range)
p_survival = (1 + np.exp(-(beta_0 + beta_1 * X + beta_2 * Y)))**-1

# plt.imshow(p_survival)
fig, ax = plt.subplots(1, 1, figsize=(4, 3))
ax.set_yscale('log')
ax.set_xlabel('effective channels per cell', fontsize=8)
ax.set_ylabel('shock rate [Hz]', fontsize=8)
ax.xaxis.set_tick_params(labelsize=8)
ax.yaxis.set_tick_params(labelsize=8)

survs = shock_data[shock_data['survival'] == True]
deaths = shock_data[shock_data['survival'] == False]
# plt.clabel(contour, inline=1, fontsize=8)

contour = ax.contourf(X, Y, p_survival, cmap='bone_r', alpha=0.5)
_ = ax.plot(survs['chan_per_cell'], survs['flow_rate'], '.', color='darkblue', ms=1,
            alpha=0.5, label='survival')
_ = ax.plot(deaths['chan_per_cell'], deaths['flow_rate'], '.', color='tomato', ms=1,
            alpha=0.5, label='death')
cbar = fig.colorbar(contour, ax=ax, ticks=[
                    0.00, 0.15, 0.3, 0.45, 0.6, 0.75, 0.90, 1.05])

cbar.ax.set_yticklabels(['0.00', '0.15', '0.30', '0.45', '0.60', '0.75', '0.90', '1'],
                        fontsize=8)
cbar.set_label('survival probability', fontsize=8)
leg = ax.legend(fontsize=8, loc='lower right')
leg.legendHandles[0]._legmarker.set_markersize(8)
leg.legendHandles[1]._legmarker.set_markersize(8)
plt.tight_layout()
plt.savefig('figS{}.pdf'.format(FIG_NO), bbox_inches='tight')
plt.savefig('figS{}.png'.format(FIG_NO), bbox_inches='tight')
