library(rstan)
library(tidyverse)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

# Read in and prepare data
synth <- all_trajectories
id <- synth[[1]]
x <- synth[[2]]
y <- synth[[3]]

J <- max(id)                  # individuals
N <- nrow(synth)

data_list <- list(
  N = N,
  J = J,
  id = id,
  x = x,
  y = y
)

# Fit model
fit <- stan(
  file = "scripts/BHM/viral_bhm.stan",
  data = data_list,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  seed = 123
)

# save fit model
#saveRDS(fit, file = "throughput/viral_bhm.rds")
