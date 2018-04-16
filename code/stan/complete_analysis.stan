data {
    // Dimensional parameters
    int<lower=1> J1; // distinct number of replicates for standard candle.
    int<lower=1> J2; // distinct number of flow rates
    int<lower=1> N1; // Total number of standard candle measurements.
    int<lower=1> N2; // Total number of shock rate measurements.
    int<lower=1, upper=J1> repl[N1]; //  Vector of trial identifier for standard candle
    int<lower=1, upper=J2> shock[N2]; // Vector of trial identifier for shock rate experiments.

    // Experimental parameters.
    vector<lower=0>[N1] A;  // Cell area for standard candle
    vector<lower=0>[N1] I_A_sc; // Intensity measurements for standard candle.
    vector<lower=0>[N2] I_A_sr; // Itensity  measurements for shock rate.
    real<lower=0> ntot;  // Mean number of channels per cell
    real<lower=0>sigma_ntot; // ERror in measurement of mean channel number.

    // Measured output
    int<lower=0, upper=1> survival[N2]; // Vector of survival identifiers.
}

parameters {
    // Standard candle parameters.
    // Hyperparameters
    real<lower=0> hyper_alpha_mu;
    real<lower=0> hyper_alpha_sigma;
    real<lower=0> hyper_A_mu;
    real<lower=0> hyper_A_sigma;

    // Low-level parameters
    vector<lower=0>[J1] alpha;
    vector<lower=0>[J1] avg_A;
    vector<lower=0>[J1] avg_I_A;
    vector<lower=0>[J1] sigma_A;
    vector<lower=0>[J1] sigma_I_A;

    // Prior for average copy number.
    real<lower=0> n_sc;

    // Shock rate parameters
    vector[J2] beta_0;
    vector[J2] beta_1;
}

transformed parameters {
     vector<lower=0>[N2] n_sr; // Channel copy number for each shocked cell.
     n_sr = I_A_sr * hyper_A_mu / hyper_alpha_mu;
}

model { 
    vector[N2] mu;
    // * Standard candle priors *//

    // Set the hyperpriors.
    n_sc ~ normal(ntot, sigma_ntot);
    hyper_A_mu ~ normal(0, 10);
    hyper_A_sigma ~ normal(0, 1);
    hyper_alpha_mu ~ normal(0, 1000);
    hyper_alpha_sigma ~ normal(0, 1000);

    // Set the low level priors.
    alpha ~ normal(hyper_alpha_mu, hyper_alpha_sigma);
    avg_A ~ normal(hyper_A_mu, hyper_A_sigma);
    sigma_A ~ normal(0, 1);


    //* Logistic regression priors *//
    beta_0 ~ normal(0, 100);
    beta_1 ~ normal(0, 100);

    // Model for standard candle.
    for (i in 1:N1) {
        A[i] ~ normal(avg_A[repl[i]], sigma_A[repl[i]]);
        avg_I_A[repl[i]] ~ normal(alpha[repl[i]] * n_sc / avg_A[repl[i]], sigma_I_A);
        I_A_sc[i] ~ normal(avg_I_A[repl[i]], sigma_I_A[repl[i]]);
    }

    // Model for logistic regression.

    for (i in 1:N2) {
        mu[i] = beta_0[shock[i]] + beta_1[shock[i]] * log(n_sr[i]);
    }
    survival ~ bernoulli_logit(mu);

}