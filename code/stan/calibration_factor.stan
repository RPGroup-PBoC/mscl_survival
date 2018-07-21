data {
    int J; // Number of unique media
    int N; // total number of measurements
    int<lower=0, upper=J> media [N]; // Media identifier
    vector<lower=0>[N] area; 
    vector<lower=0>[N] mean_intensity;
    vector<lower=0>[J] channel_mu; // Average number of channels.
    vector<lower=0>[J] channel_sig; // Error in channel measurement
}

parameters {
    real<lower=0> alpha_mu[J];
    real<lower=0> sigma[J];
    real<lower=0, upper=8> area_mu[J];
    real<lower=0> area_sigma[J];
    real<lower=0, upper=1E4> n_tot[J];
}

model {
    vector[N] mu;
    // Define the priors
    n_tot ~ normal(channel_mu, channel_sig);
    alpha_mu ~ uniform(0, 2^16);
    sigma ~ normal(0, 1);
    area_mu ~ uniform(0, 10);
    area_sigma ~ normal(0, 1);

    for (i in 1:N) {
        area[i] ~ normal(area_mu[media[i]], area_sigma[media[i]]);
        mean_intensity[i] ~ normal(alpha_mu[media[i]] * n_tot[media[i]] / area_mu[media[i]], sigma[media[i]]);
    }
    
}
