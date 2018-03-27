import numpy as np
import os
import pickle
import pystan
import pandas as pd


def load_StanModel(model_name, save_compiled=True):
    """
    Checks for a compiled Stan model given a model name and saves one if it
    does not exist.

    Parameters
    ----------
    model_name : str
        Name of the Stan model. This should end with `.stan`.
    save_compiled: bool
        If True and the compiled model does *not* exist, a compiled version
        will be saved as a `.pkl` file.

    Returns
    -------
    stan_model : pystan.StanModel instance
        The compiled stan model.

    Notes
    -----
    This function looks for all Stan models in the root `code` folder
    of this repository.
    """

    if len(model_name.split('/')[-1].split('.')) != 2:
        raise ValueError("The model name must end in .stan or .stn.")

    # Look for the compiled form of the model.
    compiled_name = model_name.split('/')[-1].split('.')[0]
    if os.path.exists('../{}.pkl'.format(compiled_name)) == True:
        print('Loading compiled model...')
        with open('../{}.pkl'.format(compiled_name), 'rb') as f:
            sm = pickle.load(f)
        print('Model loaded!')
    else:
        print('Could not find compiled code! Compiling now...')
        sm = pystan.StanModel(model_name)
        print('Compilation complete!')
        if save_compiled == True:
            print('Saving compiled model...')
            with open('../{}.pkl'.format(compiled_name), 'wb') as f:
                pickle.dump(sm, f)
            print('Model saved!')
    return sm


def chains_to_dataframe(fit, var_names=None):
    """
    Converts the generated traces from MCMC sampling to a tidy
    pandas DataFrame.

    Parameters
    ----------
    fit : pystan sampling output
        The raw MCMC output.
    var_names : list of str
        Names of desired parameters. If `None`, all parameters will be
        returned.

    Returns
    -------
    df : pandas DataFrame
        Pandas DataFrame containing all samples from the MCMC.
    """

    data_out = fit.extract()
    if var_names is None:
        _v = data_out.keys()
        var_names = []
        for v in _v:
            if v != 'lp__':
                var_names.append(v)

    for v in var_names:
        if v not in data_out.keys():
            raise ValueError("Parameter `{}` not found in index.".format(v))

    df = pd.DataFrame([])
    for k in var_names:
        shape = np.shape(data_out[k])
        if len(np.shape(data_out[k])) == 1:
            df.insert(0, k, data_out[k])
        else:
            for n in range(shape[1]):
                df.insert(0, '{}__{}'.format(k, n), data_out[k][:, n])

    logp = []
    valid_values = [v for v in data_out.keys() if v != 'lp__']
    for i in range(len(df)):
        sample = []
        for v in valid_values:
            if type(data_out[v][i]) == np.float64:
                sample.append(data_out[v][i])
            else:
                for _, k in enumerate(data_out[v][i]):
                    sample.append(k)
        logp.append(fit.log_prob(sample))
    df.insert(0, 'logp', logp)
    return df


def compute_statistics(df, var_names=None, logprob_name='logp'):
    """
    Computes the mode, hpd_min, and hpd_max from a pandas DataFrame. The value
    of the log posterior must be included in the DataFrame.

    Parameters
    ----------
    df : pandas DataFrame
        Dataframe containing MCMC samples and log posterior.
    var_names : list of str
        Name of variables desired. If None, all parameters will be
        returned.
    logprob_name : str
        Name of the log posterior column.

    Returns
    -------
    stat_df : pandas DataFrame
        DataFrame containing the mode, hpd_min, and hpd_max for each
        parameter in var_names.
    """

    # Get the vars we care about.
    if var_names == None:
        var_names = [v for v in df.keys() if v is not logprob_name]

    # Find the max of the log posterior.
    ind = np.argmax(df[logprob_name])
    if type(ind) is not int:
        ind = int(ind[0])

    # Instantiate the dataframe for the parameters.
    stat_df = pd.DataFrame([], columns=['parameter', 'mode', 'hpd_min',
                                        'hpd_max'])
    for v in var_names:
        mode = df.iloc[ind][v]
        hpd_min, hpd_max = compute_hpd(
            df[v].values.astype(float), mass_frac=0.95)
        stat_dict = dict(parameter=v, mode=mode, hpd_min=hpd_min,
                         hpd_max=hpd_max)
        stat_df = stat_df.append(stat_dict, ignore_index=True)

    return stat_df


def compute_hpd(trace, mass_frac):
    """
    Returns highest probability density region given by
    a set of samples.

    Parameters
    ----------
    trace : array
        1D array of MCMC samples for a single variable
    mass_frac : float with 0 < mass_frac <= 1
        The fraction of the probability to be included in
        the HPD.  For hreple, `massfrac` = 0.95 gives a
        95% HPD.

    Returns
    -------
    output : array, shape (2,)
        The bounds of the HPD

    Notes
    -----
    We thank Justin Bois (BBE, Caltech) for developing this function.
    http://bebi103.caltech.edu/2015/tutorials/l06_credible_regions.html
    """
    # Get sorted list
    d = np.sort(np.copy(trace))

    # Number of total samples taken
    n = len(trace)

    # Get number of samples that should be included in HPD
    n_samples = np.floor(mass_frac * n).astype(int)

    # Get width (in units of data) of all intervals with n_samples samples
    int_width = d[n_samples:] - d[:n - n_samples]

    # Pick out minimal interval
    min_int = np.argmin(int_width)

    # Return interval
    return np.array([d[min_int], d[min_int + n_samples]])
