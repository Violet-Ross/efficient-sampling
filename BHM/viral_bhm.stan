data {
  int<lower=1> N;                 // total observations
  int<lower=1> J;                 // number of individuals
  int<lower=1,upper=J> id[N];     // individual index
  vector[N] x;                    // time
  vector[N] y;                    // viral load
}

parameters {
  // Population means
  real mu_sigma;
  real mu_beta1;
  real mu_beta2;
  real mu_psi;

  // Population SDs
  real<lower=0> tau_sigma;
  real<lower=0> tau_beta1;
  real<lower=0> tau_beta2;
  real<lower=0> tau_psi;

  // Individual-level (non-centered)
  vector[J] sigma_raw;
  vector[J] beta1_raw;
  vector[J] beta2_raw;
  vector[J] psi_raw;

  // Observation noise
  real<lower=0> sigma_eps;
}

transformed parameters {
  vector[J] sigma;
  vector[J] beta1;
  vector[J] beta2;
  vector[J] psi;

  sigma = mu_sigma + tau_sigma * sigma_raw;
  beta1 = mu_beta1 + tau_beta1 * beta1_raw;
  beta2 = mu_beta2 + tau_beta2 * beta2_raw;
  psi   = mu_psi   + tau_psi   * psi_raw;
}

model {

  // Hyperpriors
  mu_sigma ~ normal(3, 1);
  mu_beta1 ~ normal(2, 1);
  mu_beta2 ~ normal(-4, 1);
  mu_psi   ~ normal(5, 1);

  tau_sigma ~ normal(0, 1);
  tau_beta1 ~ normal(0, 1);
  tau_beta2 ~ normal(0, 1);
  tau_psi   ~ normal(0, 1);

  sigma_eps ~ normal(0, 5);

  // Non-centered priors
  sigma_raw ~ normal(0, 1);
  beta1_raw ~ normal(0, 1);
  beta2_raw ~ normal(0, 1);
  psi_raw   ~ normal(0, 1);

  // Likelihood
  for (n in 1:N) {
    real mu_n;

  mu_n = sigma[id[n]]
       + beta1[id[n]] * x[n]
       + beta2[id[n]] * (x[n] - psi[id[n]]) * (x[n] > psi[id[n]]);


    y[n] ~ normal(mu_n, sigma_eps);
  }
}

