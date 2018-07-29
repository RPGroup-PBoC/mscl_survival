# -*- coding: utf-8 -*-
# %%
import sys
import numpy as np
import pandas as pd
import glob
import fcsparser
import matplotlib.pyplot as plt
import scipy.optimize
import statsmodels.tools.numdiff  as smnd
sys.path.insert(0, '../../')
import mscl.stats
import mscl.plotting
from datetime import datetime
colors = mscl.plotting.set_plotting_style()

# Include functions for gating.

# #################
def fit_2D_gaussian(df, x_val='FSC-A', y_val='SSC-A', log=False):
    '''
    This function hacks astroML fit_bivariate_normal to return the mean
    and covariance matrix when fitting a 2D gaussian fuction to the data
    contained in the x_vall and y_val columns of the DataFrame df.
    Parameters
    ----------
    df : DataFrame.
        dataframe containing the data from which to fit the distribution
    x_val, y_val : str.
        name of the dataframe columns to be used in the function
    log : bool.
        indicate if the log of the data should be use for the fit or not
    Returns
    -------
    mu : tuple.
        (x, y) location of the best-fit bivariate normal
    cov : 2 x 2 array
        covariance matrix.
        cov[0, 0] = variance of the x_val column
        cov[1, 1] = variance of the y_val column
        cov[0, 1] = cov[1, 0] = covariance of the data
    '''
    if log:
        x = np.log10(df[x_val])
        y = np.log10(df[y_val])
    else:
        x = df[x_val]
        y = df[y_val]

    # Fit the 2D Gaussian distribution using atroML function
    mu, sigma_1, sigma_2, alpha = fit_bivariate_normal(x, y, robust=True)

    # compute covariance matrix from the standar deviations and the angle
    # that the fit_bivariate_normal function returns
    sigma_xx = ((sigma_1 * np.cos(alpha)) ** 2 +
                (sigma_2 * np.sin(alpha)) ** 2)
    sigma_yy = ((sigma_1 * np.sin(alpha)) ** 2 +
                (sigma_2 * np.cos(alpha)) ** 2)
    sigma_xy = (sigma_1 ** 2 - sigma_2 ** 2) * np.sin(alpha) * np.cos(alpha)

    # put elements of the covariance matrix into an actual matrix
    cov = np.array([[sigma_xx, sigma_xy], [sigma_xy, sigma_yy]])

    return mu, cov


# #################
def gauss_interval(df, mu, cov, x_val='FSC-A', y_val='SSC-A', log=False):
    '''
    Computes the of the statistic
    (x - µx)'sum(x - µx)
    for each of the elements in df columns x_val and y_val.
    Parameters
    ----------
    df : DataFrame.
        dataframe containing the data from which to fit the distribution
    mu : array-like.
        (x, y) location of bivariate normal
    cov : 2 x 2 array
        covariance matrix
    x_val, y_val : str.
        name of the dataframe columns to be used in the function
    log : bool.
        indicate if the log of the data should be use for the fit or not.
    Returns
    -------
    statistic_gauss : array-like.
        array containing the result of the linear algebra operation:
        (x - µx)'sum(x - µx)
    '''

    # Determine that the covariance matrix is not singular
    det = np.linalg.det(cov)
    if det == 0:
        raise NameError("The covariance matrix can't be singular")

    # Compute the vector x defined as [[x - mu_x], [y - mu_y]]
    if log is True:
        x_vect = np.log10(np.array(df[[x_val, y_val]]))
    else:
        x_vect = np.array(df[[x_val, y_val]])

    x_vect[:, 0] = x_vect[:, 0] - mu[0]
    x_vect[:, 1] = x_vect[:, 1] - mu[1]

    # compute the inverse of the covariance matrix
    inv_sigma = np.linalg.inv(cov)

    # compute the operation
    interval_array = np.zeros(len(df))
    for i, x in enumerate(x_vect):
        interval_array[i] = np.dot(np.dot(x, inv_sigma), x.T)

    return interval_array


# #################
def auto_gauss_gate(df, alpha, x_val='FSC-A', y_val='SSC-A', log=False,
                    verbose=False):
    '''
    Function that applies an "unsupervised bivariate Gaussian gate" to the data
    over the channels x_val and y_val.
    Parameters
    ----------
    df : DataFrame.
        dataframe containing the data from which to fit the distribution
    alpha : float. [0, 1]
        fraction of data aimed to keep. Used to compute the chi^2 quantile
        function
    x_val, y_val : str.
        name of the dataframe columns to be used in the function
    log : bool.
        indicate if the log of the data should be use for the fit or not
    verbose : bool.
        indicate if the percentage of data kept should be print
    Returns
    -------
    df_thresh : DataFrame
        Pandas data frame to which the automatic gate was applied.
    '''
    data = df[[x_val, y_val]]
    # Fit the bivariate Gaussian distribution
    mu, cov = fit_2D_gaussian(data, log=log)

    # Compute the statistic for each of the pair of log scattering data
    interval_array = gauss_interval(data, mu, cov, log=log)

    # Find which data points fall inside the interval
    idx = interval_array <= scipy.stats.chi2.ppf(alpha, 2)

    # print the percentage of data kept
    if verbose:
        print('''
        with parameter alpha={0:0.2f}, percentage of data kept = {1:0.2f}
        '''.format(alpha, np.sum(idx) / len(df)))

    return df[idx]



#from scipy.special import erfinv
#sigmaG_factor = 1. / (2 * np.sqrt(2) * erfinv(0.5))
sigmaG_factor = 0.74130110925280102


def mean_sigma(a, axis=None, dtype=None, ddof=0, keepdims=False):
    """
    Compute mean and standard deviation for an array
    Parameters
    ----------
    a : array_like
        Array containing numbers whose mean is desired. If `a` is not an
        array, a conversion is attempted.
    axis : int, optional
        Axis along which the means are computed. The default is to compute
        the mean of the flattened array.
    dtype : dtype, optional
        Type to use in computing the standard deviation. For arrays of
        integer type the default is float64, for arrays of float types it is
        the same as the array type.
    keepdims : bool, optional
        If this is set to True, the axes which are reduced are left
        in the result as dimensions with size one. With this option,
        the result will broadcast correctly against the original `arr`.
    Returns
    -------
    mu : ndarray, see dtype parameter above
        array containing the mean values
    sigma : ndarray, see dtype parameter above.
        array containing the standard deviation
    See Also
    --------
    median_sigmaG : robust rank-based version of this calculation.
    Notes
    -----
    This routine simply calls ``np.mean`` and ``np.std``, passing the
    keyword arguments to them.  It is provided for ease of comparison
    with the function median_sigmaG()
    """
    mu = np.mean(a, axis=axis, dtype=dtype)
    sigma = np.std(a, axis=axis, dtype=dtype, ddof=ddof)

    if keepdims:
        if axis is None:
            newshape = a.ndim * (1,)
        else:
            newshape = np.asarray(a.shape)
            newshape[axis] = 1

        mu = mu.reshape(newshape)
        sigma = sigma.reshape(newshape)

    return mu, sigma


def median_sigmaG(a, axis=None, overwrite_input=False, keepdims=False):
    """
    Compute median and rank-based estimate of the standard deviation
    Parameters
    ----------
    a : array_like
        Array containing numbers whose mean is desired. If `a` is not an
        array, a conversion is attempted.
    axis : int, optional
        Axis along which the means are computed. The default is to compute
        the mean of the flattened array.
    overwrite_input : bool, optional
       If True, then allow use of memory of input array `a` for
       calculations. The input array will be modified by the call to
       median. This will save memory when you do not need to preserve
       the contents of the input array. Treat the input as undefined,
       but it will probably be fully or partially sorted.
       Default is False. Note that, if `overwrite_input` is True and the
       input is not already an array, an error will be raised.
    keepdims : bool, optional
        If this is set to True, the axes which are reduced are left
        in the result as dimensions with size one. With this option,
        the result will broadcast correctly against the original `arr`.
    Returns
    -------
    median : ndarray, see dtype parameter above
        array containing the median values
    sigmaG : ndarray, see dtype parameter above.
        array containing the robust estimator of the standard deviation
    See Also
    --------
    mean_sigma : non-robust version of this calculation
    sigmaG : robust rank-based estimate of standard deviation
    Notes
    -----
    This routine uses a single call to ``np.nanpercentile`` to find the
    quartiles along the given axis, and uses these to compute the
    median and sigmaG:
    median = q50
    sigmaG = (q75 - q25) * 0.7413
    where 0.7413 ~ 1 / (2 sqrt(2) erf^-1(0.5))
    """
    q25, median, q75 = np.nanpercentile(a, [25, 50, 75],
                                        axis=axis,
                                        overwrite_input=overwrite_input)
    sigmaG = sigmaG_factor * (q75 - q25)

    if keepdims:
        if axis is None:
            newshape = a.ndim * (1,)
        else:
            newshape = np.asarray(a.shape)
            newshape[axis] = 1

        median = median.reshape(newshape)
        sigmaG = sigmaG.reshape(newshape)

    return median, sigmaG


def sigmaG(a, axis=None, overwrite_input=False, keepdims=False):
    """
    Compute the rank-based estimate of the standard deviation
    Parameters
    ----------
    a : array_like
        Array containing numbers whose mean is desired. If `a` is not an
        array, a conversion is attempted.
    axis : int, optional
        Axis along which the means are computed. The default is to compute
        the mean of the flattened array.
    overwrite_input : bool, optional
       If True, then allow use of memory of input array `a` for
       calculations. The input array will be modified by the call to
       median. This will save memory when you do not need to preserve
       the contents of the input array. Treat the input as undefined,
       but it will probably be fully or partially sorted.
       Default is False. Note that, if `overwrite_input` is True and the
       input is not already an array, an error will be raised.
    keepdims : bool, optional
        If this is set to True, the axes which are reduced are left
        in the result as dimensions with size one. With this option,
        the result will broadcast correctly against the original `arr`.
    Returns
    -------
    median : ndarray, see dtype parameter above
        array containing the median values
    sigmaG : ndarray, see dtype parameter above.
        array containing the robust estimator of the standard deviation
    See Also
    --------
    median_sigmaG : robust rank-based estimate of mean and standard deviation
    Notes
    -----
    This routine uses a single call to ``np.nanpercentile`` to find the
    quartiles along the given axis, and uses these to compute the
    sigmaG, a robust estimate of the standard deviation sigma:
    sigmaG = 0.7413 * (q75 - q25)
    where 0.7413 ~ 1 / (2 sqrt(2) erf^-1(0.5))
    """
    q25, q75 = np.nanpercentile(a, [25, 75],
                                axis=axis,
                                overwrite_input=overwrite_input)
    sigmaG = sigmaG_factor * (q75 - q25)

    if keepdims:
        if axis is None:
            newshape = a.ndim * (1,)
        else:
            newshape = np.asarray(a.shape)
            newshape[axis] = 1

        sigmaG = sigmaG.reshape(newshape)

    return sigmaG


def fit_bivariate_normal(x, y, robust=False):
    """
    Fit bivariate normal parameters to a 2D distribution of points
    Parameters
    ----------
    x, y : array_like
        The x, y coordinates of the points
    robust : boolean (optional, default=False)
        If True, then use rank-based statistics which are robust to outliers
        Otherwise, use mean/std statistics which are not robust
    Returns
    -------
    mu : tuple
        (x, y) location of the best-fit bivariate normal
    sigma_1, sigma_2 : float
        The best-fit gaussian widths in the uncorrelated frame
    alpha : float
        The rotation angle in radians of the uncorrelated frame
    """
    x = np.asarray(x)
    y = np.asarray(y)

    assert x.shape == y.shape

    if robust:
        # use quartiles to compute center and spread
        med_x, sigmaG_x = median_sigmaG(x)
        med_y, sigmaG_y = median_sigmaG(y)

        # define the principal variables from Shevlyakov & Smirnov (2011)
        sx = 2 * sigmaG_x
        sy = 2 * sigmaG_y

        u = (x / sx + y / sy) / np.sqrt(2)
        v = (x / sx - y / sy) / np.sqrt(2)

        med_u, sigmaG_u = median_sigmaG(u)
        med_v, sigmaG_v = median_sigmaG(v)

        r_xy = ((sigmaG_u ** 2 - sigmaG_v ** 2) /
                (sigmaG_u ** 2 + sigmaG_v ** 2))

        # rename estimators
        mu_x, mu_y = med_x, med_y
        sigma_x, sigma_y = sigmaG_x, sigmaG_y
    else:
        mu_x = np.mean(x)
        sigma_x = np.std(x)

        mu_y = np.mean(y)
        sigma_y = np.std(y)

        r_xy = stats.pearsonr(x, y)[0]

    # We need to use the full (-180, 180) version of arctan: this is
    # np.arctan2(x, y) = np.arctan(x / y), modulo 180 degrees
    sigma_xy = r_xy * sigma_x * sigma_y
    alpha = 0.5 * np.arctan2(2 * sigma_xy, sigma_x ** 2 - sigma_y ** 2)

    sigma1 = np.sqrt((0.5 * (sigma_x ** 2 + sigma_y ** 2) +
                     np.sqrt(0.25 * (sigma_x ** 2 - sigma_y ** 2) ** 2 +
                     sigma_xy ** 2)))
    sigma2 = np.sqrt((0.5 * (sigma_x ** 2 + sigma_y ** 2) -
                     np.sqrt(0.25 * (sigma_x ** 2 - sigma_y ** 2) ** 2 +
                     sigma_xy ** 2)))

    return [mu_x, mu_y], sigma1, sigma2, alpha


# Start difference
T_START = 138# in sec 

# Define the data directory
files = np.sort(glob.glob('../../data/flow/RP*_.*.fcs'))
means = []
time = []
for f in files:
    m, data = fcsparser.parse(f)
    gated = auto_gauss_gate(data, 0.4, log=True)
    mean_FITC = gated['FITC-A'].mean()
    time.append(m['$ETIM']) 
    means.append(mean_FITC)


#%% Compute the proper time vector.
pattern = "%H:%M:%S"
dt = [T_START/60] 
time_init = datetime.strptime(time[0], pattern)
for i in range(1, len(time)):
    t = datetime.strptime(time[i], pattern)
    diff = t - time_init
    dt.append((T_START + diff.seconds )/ 60)

# %% Compute the growth rate
growth_data = pd.read_csv('../../data/csv/20180720_sfGFP_lb500_growth.csv')

# Trim the data to the linear region
lin_data = growth_data[(growth_data['time_min'] > 20) & (growth_data['time_min'] < 120)]
lin_data['log_A_A0'] = np.log(lin_data['od600'].values / lin_data.iloc[0]['od600'])
lin_data['adj_time'] = lin_data['time_min'] - lin_data.iloc[0]['time_min']
# Define a function to compute the grwoth rate. 
def log_post(lam, log_A_A0, time, neg=True):
    if neg == True:
        prefactor = -1
    else:
        prefactor = 1
    k = len(time)
    theo = time * lam
    return prefactor * -0.5 * np.log(k) * np.log(np.sum((log_A_A0 - theo)**2))

popt = scipy.optimize.minimize(log_post, 0.001, args=(lin_data['log_A_A0'], lin_data['adj_time']), 
                                method='powell')
hess = smnd.approx_hess(np.array([float(popt.x)]), log_post, args=(lin_data['log_A_A0'], lin_data['adj_time'], False))
cov = -np.linalg.inv(hess)
std = np.sqrt(cov)
mean_lam = float(popt.x)
err = np.sqrt(cov)[0][0]


# %%

# Instantiate the figure
fig, ax = plt.subplots(1, 2, figsize=(6, 3))

# Format axes 
for a in ax:
    a.xaxis.set_tick_params(labelsize=8)
    a.yaxis.set_tick_params(labelsize=8)
ax[0].set_ylim([0.3, 1.15])
ax[1].set_xlim([-5, 80])

# Add labels
ax[0].set_xlabel('time [min]', fontsize=8)
ax[0].set_ylabel('normalized fluorescence', fontsize=8)
ax[1].set_xlabel('time [min]', fontsize=8)
ax[1].set_ylabel('$\log (A / A_0)$')
ax[0].text(-0.2, 1.0, '(A)', fontsize=8, transform=ax[0].transAxes)
ax[1].text(-0.24, 1.0, '(B)', fontsize=8, transform=ax[1].transAxes)
ax[0].set_xlabel('time [min]')
ax[0].set_ylabel('normalized fluorescence')

# Plot the maturation data
_ = ax[0].plot(dt[:-1], np.array(means)/means[-1], '-o', lw=1, ms=3, color=colors['red'], label='measured\nintensity')
_ = ax[0].fill_betweenx(np.linspace(0, 1.15, 500), 30, 40, color=colors['pale_yellow'], zorder=-1, 
                label='observed\n growth rate')

time_range = np.linspace(-5, 80, 500)
_ = ax[1].plot(time_range, mean_lam * time_range, 'k-', lw=0.5, label='best fit')
_ = ax[1].fill_between(time_range, (mean_lam - std)[0][0] * time_range,
                    (mean_lam + std)[0][0] * time_range, color='slategray',
                    alpha=0.5, label='__nolegend__')
_ = ax[1].plot(lin_data['adj_time'], lin_data['log_A_A0'], 'o', color=colors['red'], ms=4,
               label='data')
leg = ax[1].legend(loc='upper left', fontsize=8, title='$t_{double} = 35 \pm 1$ min') 
leg.get_title().set_fontsize(8)
plt.tight_layout()

plt.savefig('../../figs/figRX_maturation_time.pdf')
plt.savefig('../../figs/figRX_maturation_time.png', dpi=300)

