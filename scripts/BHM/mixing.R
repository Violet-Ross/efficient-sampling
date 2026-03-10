library(rstan)
library(tidyverse)
library(bayesplot)
library(patchwork)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

fit <- readRDS("throughput/viral_bhm.rds")

posterior_array <- rstan::extract(
  fit,
  permuted = FALSE
)

p_trace <- mcmc_trace(
  posterior_array,
  pars = "mu_alpha"
)

p_acf <- mcmc_acf(
  posterior_array,
  pars = "mu_alpha",
  lags = 30
)

mixing_plot <- p_trace + p_acf

ggsave("figures/mixing.png", plot = mixing_plot)
