library(rstan)
library(tidyverse)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

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

fit <- stan(
  file = "viral_bhm.stan",
  data = data_list,
  chains = 4,
  iter = 2000,
  warmup = 1000,
  seed = 123
)

saveRDS(fit, file = "stan_fit.rds")

print(fit, pars = c("mu_alpha", "mu_beta1",
                    "mu_beta2", "mu_psi",
                    "sigma_eps"))


# ---- Population trajectory evaluated at posterior means ----

# Extract posterior draws
post <- rstan::extract(fit)

# Time grid
x_grid <- seq(min(x), max(x), length.out = 200)
n_draws <- length(post$mu_alpha)

# Compute trajectory for all draws
curve_draws <- sapply(1:n_draws, function(i) {
  post$mu_alpha[i] +
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

