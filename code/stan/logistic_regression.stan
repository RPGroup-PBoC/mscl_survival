data { 
    int<lower=0> J; // Number of distinct data sets
    int<lower=0> N; // Number of data points
    int<lower=0, upper=J> trials[N]; // Sequence of identifiers for each measurement
    vector<lower=0>[N] predictor; // Vector for predictor value
    vector<lower=0>[N] predictor_err; // Statistical error for predictor
    int<lower=0, upper=1> output[N]; // Boolean output vector
}

parameters {
    real beta_0[J]; // Intercept for each trial
    real beta_1[J]; // Slope for each trial
    vector<lower=0>[N] predictor_mu; // The most-likely value for the predictor
}

model {
    // Priors for slope and intercept. Assign weakly informative distributions
    beta_0 ~ normal(0, 100);
    beta_1 ~ normal(0, 100);
    // Compute the most likely value for the predictor value
    predictor_mu ~ normal(predictor, predictor_err);
    // Loop through each trial and compute the likelihood using appropriate params.
    for (i in 1:N) {
        output ~ bernoulli_logit(beta_0[trials[i]] + beta_1[trials[i]] * predictor_mu);
    }
}