data {
    // Dimensional parameters
    int<lower=1> J; // Number of replicates
    int<lower=1> N; // Number of individual measurements
    int<lower=1, upper=J> rep[N]; // vector of replicate idenfiiers

    // Experimental parameters
    vector<lower=0>[N] A; // Area
    vector<lower=0>[N] I_A; // Areal intensity
    real<lower=0> ntot; // mean number of channels per cell    
    real<lower=0> sigma_ntot; // Error in measurement of the mean
}

parameters { 

    // Define the low level parameters
    real<lower=0> alpha[J]; // Calibration factor.
    real<lower=0> avg_A[J]; // Average area 
    real<lower=0> avg_I_A[J]; // Average areal intensity
    real<lower=1E-6> sigma_A[J]; // Error in area
    real<lower=1E-6> sigma_I_A[J]; // Error in areal intensity

    // Define the hyper parameters for the means.
    real<lower=0> hyper_alpha_mu;
    real<lower=0> hyper_A_mu;
    real<lower=1E-6> hyper_alpha_mu_sigma;
    real<lower=1E-6> hyper_A_mu_sigma;
    real<lower=1E-6> hyper_I_A_mu_sigma;

    // Define the hyper parameters for the variance.
    real<lower=0> hyper_alpha_sigma;
    real<lower=0> hyper_A_sigma;
    real<lower=0> hyper_I_A_sigma;
    real<lower=1E-6> hyper_alpha_sigma_sigma;
    real<lower=1E-6> hyper_A_sigma_sigma;
    real<lower=1E-6> hyper_I_A_sigma_sigma;

    // Define the prior for the average copy number.
    real <lower=0> n;
}

model {

    n ~ normal(ntot, sigma_ntot);

    // Set the hyper priors for the means.
    hyper_A_mu ~ normal(0, 1000);
    hyper_alpha_mu_sigma ~ normal(0, 100);
    hyper_A_mu_sigma ~ normal(0, 100);
    hyper_I_A_mu_sigma ~ normal(0, 100);

    // Set the hyperpriors for the variance.
    hyper_A_sigma ~ normal(0, 100);
    hyper_I_A_sigma ~ normal(0, 100);
    hyper_A_sigma_sigma ~ normal(0, 100);
    hyper_I_A_sigma_sigma ~ normal(0, 100);

    // Define the low level priors.
    alpha ~ normal(hyper_alpha_mu, hyper_alpha_mu_sigma);
    avg_A ~ normal(hyper_A_mu, hyper_A_mu_sigma);
    sigma_A ~ normal(hyper_A_sigma, hyper_A_sigma_sigma);
    sigma_I_A ~ normal(hyper_I_A_sigma, hyper_I_A_sigma_sigma);

     
    for (i in 1:N) {
        A[i] ~ normal(avg_A[rep[i]], sigma_A[rep[i]]);
        avg_I_A[rep[i]] ~ normal(alpha[rep[i]]* n / avg_A[rep[i]] , sigma_I_A[rep[i]]);
        I_A[i] ~ normal(avg_I_A[rep[i]], sigma_I_A[rep[i]]);
    }
    
}