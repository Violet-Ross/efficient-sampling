library(rstan)
library(tidyverse)
library(bayesplot)
library(patchwork)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Read in and prepare data
synth <- read_csv("BHM/synthetic.csv")
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
  iter = 1000,
  warmup = 500,
  seed = 123
)


# --- 2. Extract posterior draws as array ---
posterior_array <- rstan::extract(
  fit,
  permuted = FALSE
)

# --- 3. Trace plot (mixing) ---
p_trace <- mcmc_trace(
  posterior_array,
  pars = "mu_alpha"
)

# --- 4. ACF plot ---
p_acf <- mcmc_acf(
  posterior_array,
  pars = "mu_alpha",
  lags = 30
)

# --- 5. Combine side by side ---
mixing_plot <- p_trace + p_acf

ggsave("figures/mixing.png", plot = mixing_plot)
