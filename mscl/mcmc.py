import pymc3 as pm
import theano.tensor as tt
import numpy as np
import pandas as pd


def _log_prior_trace(trace, model):
    """
    Computes the contribution of the log prior to the log posterior.

    Parameters
    ----------
    trace : PyMC3 trace object.
        Trace from the PyMC3 sampling.
    model : PyMC3 model object
        Model under which the sampling was performed

    Returns
    -------
    log_prior_vals : nd-array
        Array of log-prior values computed elementwise for each point in the
        trace.

    Notes
    -----
    This function was modified from one produced by Justin Bois.
    http://bebi103.caltech.edu
    """
    # Iterate through each trace.
    try:
        points = trace.points()
    except:
        points = trace

    # Get the unobserved variables.
    priors = [var.logp for var in model.unobserved_RVs if type(
        var) == pymc3.model.FreeRV]

    def logp_vals(pt):
        if len(model.unobserved_RVs) == 0:
            return pm.theanof.floatX(np.array([]), dtype='d')

        return np.array([logp(pt) for logp in priors])

    # Compute the logp for each value of the prior.
    log_prior = (logp_vals(pt) for pt in points)
    return np.stack(log_prior)


def _log_post_trace(trace, model):
    R"""
    Computes the log posterior of a PyMC3 sampling trace.

    Parameters
    ----------
    trace : PyMC3 trace object
        Trace from MCMC sampling
    model: PyMC3 model object
        Model under which the sampling was performed.

    Returns
    -------
    log_post : nd-array
        Array of log posterior values computed elementwise for each point in
        the trace

    Notes
    -----
    This function was modified from one produced by Justin Bois
    http://bebi103.caltech.edu
    """

    # Compute the log likelihood. Note this is improperly named in PyMC3.
    log_like = pm.stats._log_post_trace(trace, model).sum(axis=1)

    # Compute the log prior
    log_prior = _log_prior_trace(trace, model)

    return (log_piror.sum(axis=1) + log_like)


def trace_to_dataframe(trace, model):
    R"""
    Converts a PyMC3 sampling trace object to a pandas DataFrame

    Parameters
    ----------
    trace, model: PyMC3 sampling objects.
        The MCMC sampling trace and the model context.

    Returns
    -------
    df : pandas DataFrame
        A tidy data frame containing the sampling trace for each variable  and
        the computed log posterior at each point.
    """

    # Use the Pymc3 utilitity.
    df = pm.trace_to_dataframe(trace)

    # Include the log prop
    df['logp'] = _log_post_trace(trace, model)
    return df


def trace_to_dataframe(trace, model):
    """
    Converts a PyMC3 sampling trace object to a pandas DataFrame
    """

    # Use the Pymc3 utilitity.
    df = pm.trace_to_dataframe(trace)

    # Include the log prop
    df['logp'] = pm.stats._log_post_trace(trace, model).sum(axis=1)
    return df


def compute_statistics(df, varnames=None, logprob_name='logp'):
    """
    Computes the mode, hpd_min, and hpd_max from a pandas DataFrame. The value
    of the log posterior must be included in the DataFrame.
    """

    # Get the vars we care about.
    if varnames is None:
        varnames = [v for v in df.keys() if v is not 'logp']

    # Find the max of the log posterior.
    ind = np.argmax(df[logprob_name])
    if type(ind) is not int:
        ind = ind[0]

    # Instantiate the dataframe for the parameters.
    stat_df = pd.DataFrame([], columns=['parameter', 'mode', 'hpd_min',
                                        'hpd_max'])
    for v in varnames:
        mode = df.iloc[ind][v]
        hpd_min, hpd_max = compute_hpd(df[v].values, mass_frac=0.95)
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


def logistic(val):
    """
    Computes the logistic function using Theano. Logistic function is
    defined as
        logistic = (1 + e^(-val))^-1
    """
    return (1 + tt.exp(-val))**-1


class MarginalizedNormal(pm.Continuous):
    """
    A bivariate Normal distribution after marginalization of sigma.

    g(µ, k| y) = ((y - μ)^2)^(-k/2)

    Parameters
    ----------
    mu : PyMC3 RV object
        The mean of the components of the distribution.

    """

    def __init__(self, mu=None, *args, **kwargs):
        super(MarginalizedNormal, self).__init__(*args, **kwargs)
        self.mu = mu = pm.theanof.floatX(tt.as_tensor_variable(mu))
        self.median = mu
        self.mode = mu
        self.mean = mu

    def logp(self, values):
        k = values.shape[-1]
        mu = self.mu
        return -0.5 * k * tt.log(tt.sum((values - mu)**2))


class GammaApproxBinomial(pm.Continuous):
    """
    An approximation of the Binomial distribution where the Binomial
    coefficient is calculated via gamma functions,

    n! = nΓ(n) = Γ(n + 1)

    Parameters
    ----------
    n, p : PyMC3 RV objects
        The number of Bernoulli trials (n) and the probability of success
        (p). Probability must be on the range p in [0, 1]
    """

    def __init__(self, n=None, p=None, *args, **kwargs):

        # Ensure parameters are defined and properly bounded.
        if n is None or p is None:
            raise RuntimeError('vars {0} and {1} must be defined'.format(n, p))
        super(GammaApproxBinomial, self).__init__(*args, **kwargs)
        self.n = n = pm.theanof.floatX(tt.as_tensor_variable(n))
        self.p = p = pm.theanof.floatX(tt.as_tensor_variable(p))

        # Define the testvals.
        self.mean = n * p

    def logp(self, value):
        p = self.p
        n = self.n
        binomial_coeff = tt.gammaln(n + 1) - tt.gammaln(value + 1) -\
            tt.gammaln(n - value + 1)
        prob = value * tt.log(p) + (n - value) * tt.log(1 - p)
        return binomial_coeff + prob


class Jeffreys(pm.Continuous):
    """
    Jeffreys prior for a scale parameter.

    Parameters
    ----------
    lower : float, > 0
        Minimum value the variable can take.
    upper : float, > `lower`
        Maximum value the variable can take.
    Returns
    -------
    output : pymc3 distribution
        Distribution for Jeffreys prior.

    Notes
    -----
    This class was adopted from Justin Bois
    github.com/justinbois/bebi103
    """

    def __init__(self, lower=None, upper=None, transform='interval',
                 *args, **kwargs):
        # Check inputs
        if lower is None or upper is None:
            raise RuntimeError('`lower` and `upper` must be provided.')

        if transform == 'interval':
            transform = pm.distributions.transforms.interval(lower, upper)
        super(Jeffreys, self).__init__(transform=transform, *args, **kwargs)
        self.lower = lower = pm.theanof.floatX(tt.as_tensor_variable(lower))
        self.upper = upper = pm.theanof.floatX(tt.as_tensor_variable(upper))

        self.mean = (upper - lower) / tt.log(upper / lower)
        self.median = tt.sqrt(lower * upper)
        self.mode = lower

    def logp(self, value):
        lower = self.lower
        upper = self.upper
        return pm.distributions.dist_math.bound(
            -tt.log(tt.log(upper / lower)) - tt.log(value),
            value >= lower, value <= upper)


def ReparameterizedNormal(name=None, mu=None, sd=None, shape=1):
    """
    A reparameterized normal distribution.

    Parameters
    ----------
    name :  string
        The name of the RV. The reparameterized version will have this name prepended with "offset_"
    mu : float
        Mean of the normal distribution.
    sd: float
        The standard deviation if the distribtion.
    shape : int
        The shape of the RV. Default is 1
    """
    if name is None:
        raise RuntimeError("`name` must be provided.")
    if mu is None:
        raise RuntimeError("`mu` must be provided.")
    if sd is None:
        raise RuntimeError("`sd` must be provided.")
    if type(name) is not str:
        raise TypeError(
            "expected type(name) to be string, got {0}.".format(type(name)))

    # Compute the offset.
    offset_var = pm.Normal('offset_{0}'.format(name), mu=0, sd=1, shape=shape)

    # Define the reparameterized variable.
    var = pm.Deterministic(name, mu + offset_var * sd)
    return var
