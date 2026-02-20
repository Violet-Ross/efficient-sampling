library(rstan)
library(tidyverse)
library(patchwork)

fit <- readRDS("stan_fit.rds")

post <- rstan::extract(fit)
true_params <- read_csv("true_population_params.csv")

actual_alpha <- rnorm(1000000, true_params$mu_alpha, true_params$sigma_alpha)
p1 <- ggplot() +
  geom_density(aes(x = post$alpha, color = "posterior"),
               key_glyph = "path") +
  geom_density(aes(x = actual_alpha, color = "true distribution"),
               key_glyph = "path") +
  labs(color = NULL, x = "sigma")


actual_beta1 <- rnorm(1000000, true_params$mu_beta1, true_params$sigma_beta1)
p2 <- ggplot() +
  geom_density(aes(post$beta1, color = "posterior"),  key_glyph = "path") +
  geom_density(aes(actual_beta1, color = "true distribution"),  key_glyph = "path") +
  labs(color = NULL, x = "beta1")

actual_beta2 <- rnorm(1000000, true_params$mu_beta2, true_params$sigma_beta2)
p3 <- ggplot() +
  geom_density(aes(x = post$beta2, color = "posterior"),
               key_glyph = "path") +
  geom_density(aes(x = actual_beta2, color = "true distribution"),
               key_glyph = "path") +
  labs(color = NULL, x = "beta2")

actual_psi <- rnorm(1000000, true_params$mu_psi, true_params$sigma_psi)
p4 <- ggplot() +
  geom_density(aes(x = post$psi, color = "posterior"),
               key_glyph = "path") +
  geom_density(aes(x = actual_psi, color = "true distribution"),
               key_glyph = "path") +
  labs(color = NULL, x = "psi")

plot_500_individs <- 
  (p1 | p2) / (p3 | p4) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("plot_500_individs.png", plot = plot_500_individs)

plot_500_individs


