# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
import statsmodels.tools.numdiff as smnd
import scipy.optimize
import glob


def ecdf(data):
    return np.sort(data), np.arange(0, len(data)) / len(data)


def density_binning(data, groupby='shock_group', bin_width=10,
                    min_events=20, input_key='channel_density',
                    output_key='survival'):
    """
    Bins an array of data by a given density.

    Parameters
    ----------
    data : pandas DataFrame
        DataFrame containing data with computed channel density,
        survival classifier, and shock speed designation.
    groupby : list of strings.
        Keys by which to group the survival data. Default is 'shock_group'
    channel_bin : float or int
        Bin width for channel density. Default is 10 channels per unit area.
    min_cells : int
        Minimum number of cells to consider for each bin.
    channel_key : string
        Column name for channel density. Default is 'channel_density'.
    survival_key : string
        Column name for survival identifier. Default is 'survival'.

    Returns:
    --------
    bin_data : pandas DataFrame
        Data frame with binned data.
    """

    # Set the bounds for the bins.
    lower_bound = 0
    upper_bound = int(data[input_key].max())
    bins = np.arange(lower_bound, upper_bound + bin_width, bin_width)

    # Sort the data by channel density
    sorted_data = data.sort_values(by=input_key)

    # Partition into the bins.
    sorted_data['bin_number'] = 0
    sorted_data.reset_index(inplace=True)
    for i in range(1, len(bins) - 1):
        # Assign bin numbers based on channel density
        inds = (sorted_data[input_key] >= bins[i - 1]
                ) & (sorted_data[input_key] < bins[i + 1])
        sorted_data.loc[inds,  'bin_number'] = i

    # Ensure that the bin numbering scheme is sequential.
    bin_data = sorted_data.copy()
    grouped = bin_data.groupby(groupby)
    for g, d in grouped:
        seq_change = {}
        bin_nos = d['bin_number'].unique()
        for i, b in enumerate(bin_nos):
            bin_data.loc[(bin_data['bin_number'] == b) &
                         (bin_data[groupby] == g), 'bin_number'] = i

    # Regroup the data and congeal bins to a minimum cell number.
    grouped = bin_data.groupby(groupby)
    for g, d in grouped:
        # Group by bin number.
        _grouped = d.groupby('bin_number')[output_key].count()

        # Find those less than the minimum cell number.
        bin_counts = _grouped.to_dict()
        low_bins = _grouped[_grouped < min_events].to_dict()

        # Get just the bin numbers.
        bin_nos = list(low_bins.keys())

        # Identify the edges of sequential bins with low cell counts.
        sequential = np.where(np.diff(bin_nos) > 1)[0]
        if (len(sequential) == 0) & (len(bin_nos) != 0):
            paired = [bin_nos]
        else:
            # Need to do fancy indexing here so it returns even single bins.
            paired = [bin_nos[:sequential[0] + 1]]
            _paired = ([bin_nos[sequential[j - 1] + 1:sequential[j] + 1]
                        for j in range(1, len(sequential))])
            for _p in _paired:
                paired.append(_p)
            paired.append(bin_nos[sequential[-1] + 1:])

        # Loop through each pair and determine if they can meet the minimum.
        change_bins = {}

        for i, pair in enumerate(paired):
            if len(pair) > 1:
                summed = np.sum([bin_counts[p] for p in pair])
                if summed >= min_cells:
                    for z in pair:
                        change_bins[z] = pair[0]
                else:
                    for z in pair:
                        change_bins[z] = pair[0] - 1
            else:
                # Deal with edge cases of first and last bin.
                if pair[0] == 1:
                    change_bins[pair[0]] = pair[0] + 1
                elif pair[0] == sorted_data['bin_number'].max():
                    change_bins[pair[0]] = pair[0] - 1

        # Loop through the changed bins and change the value of the bin
        # number in the original dataframe.
        keys = change_bins.keys()
        for key in keys:
            bin_data.loc[(bin_data[groupby] == g) &
                         (bin_data['bin_number'] == key),
                         'bin_number'] = change_bins[key]
    return bin_data


def compute_survival_stats(df, groupby=['shock_group', 'bin_number']):
    """
    Computes the statistics of survival probabilitiy, number of cells, and
    binomial error given a dataframe with binned events. This should be used
    as an apply function on a pandas groupby method.
    """
    def binomial_probability(df):
        n = np.sum(df == True)
        N = len(df)
        return n / N

    def binomial_err(df):
        n = np.sum(df == True)
        N = len(df)
        return np.sqrt(n * (N - n) / N**3)

    def _compute_stats(df):
        stats_dict = dict(prob=binomial_probability(df['survival']),
                          err=binomial_err(df['survival']),
                          mean_chan=df['channel_density'].mean(),
                          n_cells=len(df), n_suv=np.sum(df['survival']))
        return pd.Series(stats_dict)
    grouped = df.groupby(groupby).apply(_compute_stats)
    return pd.DataFrame(grouped).reset_index()


def deterministic_log_post(alpha, I1, I2, p=0.5, neg=True):
    """
    Computes the log posterior for the deterministic measurement
    of the calibration factor α.

    Parmeters
    ---------
    alpha : float or int
        The calibration factor in units of A.U per molecule. This parameter
        must be positive.
    I1, I2 : arrays of floats or ints.
        The intensities of the two daugher cells. Each should be provided
        individually as 1d-arrays.
    p : float
        The probability of partitioning proteins into the first daughter cell.
        Defalt value is 0.5.
    neg : bool
        If True, the negative log posterior will be returned.

    Returns
    -------
    log_post: 1d-array
        The value of the log posterior at a given value of α.
    """
    # Ensure positivity of samples.
    if alpha < 0:
        return -np.inf
    if p < 0 or p > 1:
        raise RuntimeError('p must be between 0 and 1.')

    # Set the prior
    k = len(I1)
    prior = -k * np.log(alpha)

    # Determine the protein copy numbers deterministically.
    n1 = I1 / alpha
    n2 = I2 / alpha
    ntot = n1 + n2

    # Compute the binomial coefficient using log gamma functions.
    binom = scipy.special.gammaln(ntot + 1).sum() -\
        scipy.special.gammaln(n1 + 1).sum() - \
        scipy.special.gammaln(n2 + 1).sum()

    # Compute the probability factor of the binomial.
    prob = n1.sum() * np.log(p) + n2.sum() * np.log(1 - p)

    # Determine if the negative log posterior is desired.
    if neg is True:
        prefactor = -1
    else:
        prefactor = 1

    # Return the desired quantity.
    return prefactor * (prior + prob + binom)


def log_posterior_biexp(params, time, data, neg=True):
    """
    Computes the log posterior for a single exponential decay.

    Parameters
    ----------
    params : list or tuple
        The parameter values for the model. These should be in the
        order of beta, I_0_1, tau_1, I_0_2, tau_2.
    time: list or 1d-array
        The time points overwhich the bleaching occurs.
    data : list or 1d-array
        The bleaching data
    neg : bool
        If True, the negative log posterior is returned. Default is True.

    Returns
    -------
    log_post : float
        The value of the log posterior given the parameter inputs.
    """
    if neg is True:
        prefactor = -1
    else:
        prefactor = 1

    if (params < 0).any():
        return prefactor * -np.inf

    beta, I_0_1, tau_1, I_0_2, tau_2 = params

    k = len(data)
    mu = beta + I_0_1 * np.exp(-time / tau_1) + I_0_2 * np.exp(-time / tau_2)
    log_post = -np.log(tau_1) - np.log(tau_2) - (k / 2) * \
        np.log(np.sum((mu - data)**2))
    return prefactor * log_post


def estimate_calibration_factor(I1, I2, p=0.5):
    """
    Estimates the optimal value of α for a given data set by minimization.

    Parameters
    ----------
    I1, I2 : 1d-arrays
        The intensities of the two daughter cells. These should be provided
        individually as 1d-arrays.
    p : float
        The probability of paritioning into the first daughter cell. Default
        value is even, 0.5.

    Returns
    -------
    alpha_opt : float
        Best-fit value for the value of α in units of A.U. per molecule.
    sd : float
        The standard deviation of the the gaussian approximation of the
        posterior.
    """

    # Minimize the negative log posterior.
    popt = scipy.optimize.minimize_scalar(deterministic_log_post,
                                          args=(I1, I2, p, True))
    alpha_opt = popt.x

    # Compute the hessian.
    hess = smnd.approx_hess([alpha_opt], deterministic_log_post,
                            args=(I1, I2, p, False))
    cov = -np.linalg.inv(hess)
    sd = np.sqrt(cov[0])

    return [alpha_opt, sd[0]]
