import pymc3 as pm
import theano.tensor as tt


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
        mu_array = self.mu_array
        return -0.5 * k * tt.log(tt.sum((values - mu_array)**2))


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
