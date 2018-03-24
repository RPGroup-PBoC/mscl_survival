/**
* Heriarchical Losgistic Regression
* --------------------------------------------------------------------------
* This Stan model performs Bayesian parameter estimation for the coefficients
* of a single regressor logistic regression including statistical error for
* each single cell channel number.
*
* Author: Griffin Chure
* Creation Date: March 18, 2018
* License: MIT
* Contact: gchure@caltech.edu
**/

data {
    int<lower=0> N; // Total number of sample measurements.
    int<lower=0, upper=1> output[N]; // Categorial input vector.
    vector<lower=0>[N] predictor;  // Computed effective channel number.
    vector<lower=0>[N] predictor_err; // Statistical error in channel number.
    }

parameters {
    real beta_0; // Intercept - log odds of survival with no channels.
    real beta_1; // Slope - increase in log odds from +1 increase in channels.
    vector<lower=0>[N] predictor_mu; // Single-cell channel number.
    }

model {
    beta_0 ~ normal(0, 100); // Weakly informative prior for intercept
    beta_1 ~ normal(0, 100); // Weakly informative prior for slope

    // Informative prior for channel number with provided error.
    predictor_mu ~ normal(predictor, predictor_err);

    // Likelihood as a logit Bernoulli distribution.
    output ~ bernoulli_logit(beta_0 + beta_1 * predictor);
  }
