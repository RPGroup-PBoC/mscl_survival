data {
    // Dimensional parameters
    int<lower=1> J; // Number of unique classes
    int<lower=1> N; // Number of measurements
    int<lower=1, upper=J> trial[N]; // Vector of trial identifiers.

    // Experimental parameters.
    vector[N] I_A; // Vector of measured areal intensities
    real alpha; // Calibration factor
    real alpha_sigma; // Error in alpha
    real A; // Area  
    real A_sigma; // Error in area.

    // Measured output
    int survival[N]; // Survival
}

parameters {
    real beta_0[J]; // Intercept
    real beta_1[J]; // Slope
    real<lower=0> alpha_mu;
    real<lower=0> A_mu;
    }

transformed parameters {
    vector<lower=0>[N] n;
    n = I_A * A_mu / alpha_mu;
}

model {
    vector[N] mu;
    alpha_mu ~ normal(alpha, alpha_sigma);
    A_mu ~ normal(A, A_sigma);
    beta_0 ~ normal(0, 100);
    beta_1 ~ normal(0, 100);
    for (i in 1:N) {
        mu[i] = beta_0[trial[i]] + beta_1[trial[i]] * log(n[i]);
    }
    survival ~ bernoulli_logit(mu);
}
