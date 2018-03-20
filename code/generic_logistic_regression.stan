/**
* Heriarchical Losgistic Regression
* --------------------------------------------------------------------------
* This Stan model performs Bayesian parameter estimation for the coefficients
* of a single regressor logistic regression neglecting any statistical error.
*
* Author: Griffin Chure
* Creation Date: March 19, 2018
* License: MIT
* Contact: gchure@caltech.edu
**/

data {
    int<lower=0> N; // Total number of sample measurements.
    int<lower=0, upper=1> output[N]; // Categorial input vector.
    vector<lower=0>[N] predictor;  // Computed effective channel number.
    }

parameters {
    real beta_0; // Intercept - log odds of survival with no channels.
    real beta_1; // Slope - increase in log odds from +1 increase in channels.
    }

model {
    beta_0 ~ normal(0, 100); // Weakly informative prior for intercept
    beta_1 ~ normal(0, 100); // Weakly informative prior for slope

    // Likelihood as a logit Bernoulli distribution.
    output ~ bernoulli_logit(beta_0 + beta_1 * predictor);
    }
