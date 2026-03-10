library(rstan)
library(tidyverse)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Read in and prepare data
synth <- read_csv("synthetic.csv")
id <- synth[[1]]
x <- synth[[2]]
y <- synth[[3]]

J <- max(id)                  # individuals
n_per <- sum(id == 1)         # observations per individual
N <- J * n_per

data_list <- list(
  N = N,
  J = J,
  id = id,
  x = x,
  y = y
)

# Fit model
fit <- stan(
  file = "viral_bhm.stan",
  data = data_list,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  seed = 123
)

# save fit model
saveRDS(fit, file = "stan_fit.rds")

# print the posteriors
print(fit, pars = c("mu_alpha", "mu_beta1",
                    "mu_beta2", "mu_psi",
                    "sigma_eps"))