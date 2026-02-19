library(rstan)
library(tidyverse)
library(patchwork)

fit <- readRDS("stan_fit.rds")

post <- rstan::extract(fit)

actual_alpha <- rnorm(1000000, mu_alpha, sigma_alpha)
p1 <- ggplot() +
  geom_density(aes(x = post$sigma, color = "posterior"),
               key_glyph = "path") +
  geom_density(aes(x = actual_alpha, color = "true distribution"),
               key_glyph = "path") +
  labs(color = NULL, x = "sigma")


actual_beta1 <- rnorm(1000000, mu_beta1, sigma_beta1)
p2 <- ggplot() +
  geom_density(aes(post$beta1, color = "posterior"),  key_glyph = "path") +
  geom_density(aes(actual_beta1, color = "true distribution"),  key_glyph = "path") +
  labs(color = NULL, x = "beta1")

actual_beta2 <- rnorm(1000000, mu_beta2, sigma_beta2)
p3 <- ggplot() +
  geom_density(aes(x = post$beta2, color = "posterior"),
               key_glyph = "path") +
  geom_density(aes(x = actual_beta2, color = "true distribution"),
               key_glyph = "path") +
  labs(color = NULL, x = "beta2")

actual_psi <- rnorm(1000000, mu_psi, sigma_psi)
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


# ---- Bonus Plot: Population trajectory evaluated at posterior means ----

# Extract posterior draws
post <- rstan::extract(fit)

# Time grid
x_grid <- seq(min(x), max(x), length.out = 200)
n_draws <- length(post$mu_sigma)

# Compute trajectory for all draws
curve_draws <- sapply(1:n_draws, function(i) {
  post$mu_sigma[i] +
    post$mu_beta1[i] * x_grid +
    post$mu_beta2[i] * pmax(0, x_grid - post$mu_psi[i])
})

# Convert to long data frame
df <- as.data.frame(curve_draws)
df$x <- x_grid
df_long <- pivot_longer(df, cols = -x, names_to = "draw", values_to = "y")

# Compute mean and 95% CI
df_summary <- df_long %>%
  group_by(x) %>%
  summarize(mean = mean(y),
            lower = quantile(y, 0.025),
            upper = quantile(y, 0.975))

# Plot with ggplot2
ggplot(df_summary, aes(x = x, y = mean)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "blue", alpha = 0.2) +
  geom_line(color = "blue", size = 1.5) +
  #geom_point(data = synth %>% filter(PersonID == 1), aes(x = TestDateIndex, y = CtT1, group = factor(PersonID))) +
  labs(x = "Time", y = "Viral Load",
       title = "Population Viral Trajectory with 95% CI") +
  theme_minimal()


