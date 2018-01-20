import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import mscl.plotting
colors = mscl.plotting.set_plotting_style()


# Load the survival error estimate data.
err_estimate = pd.read_csv(
    '/Users/gchure/Desktop/VanDeBerg_2016_Supp_Survival.csv')
surv_data = pd.read_csv('/Users/gchure/Desktop/Survival_Plot_Poolman.csv')
chan_err = pd.read_csv('/Users/gchure/Desktop/Poolman_Err_Data.csv')

# Compute the error.
y_err = err_estimate.values[:, 1]
avg_yerr_perc = np.mean(np.diff(np.sort(y_err)) / y_err[1]) * 100

# Set up the plot.
fig, ax = plt.subplots(1, 1)
_ = ax.plot(surv_data['chan_per_cell'],
            surv_data['surv_perc'], 'o', color=colors['red'])
_ = ax.errorbar(surv_data['chan_per_cell'], surv_data['surv_perc'],
                yerr=np.ones(len(surv_data)) * avg_yerr_perc, xerr=chan_err['chan_per_cell'] - surv_data['chan_per_cell'],
                linestyle='none', color=colors['red'])

ax.set_ylim([0, 130])
