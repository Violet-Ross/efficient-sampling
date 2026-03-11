library(rstan)
library(tidyverse)

options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

true_params <- read_csv("throughput/true_population_params.csv")
synth <- read_csv("throughput/synthetic.csv")
y_true <- synth[[3]]

max_samples = 500
step_size = 50
n_values <- seq(from = step_size, to = max_samples, by = step_size)
steps <- length(n_values)

rmse_by_size <- matrix(NA, steps, 5)
coverage_by_size <- matrix(NA, steps, 5)

model <- stan_model("scripts/BHM/viral_bhm.stan")

for(i in seq_along(n_values)){
  
  n <- n_values[i]
  print(paste0("Fitting model to ", n, " samples"))
  
  sample <- data[data$ids <= n, ]
  
  id <- sample[[1]]
  x <- sample[[2]]
  y <- sample[[3]]
  
  J <- length(unique(id))
  n_per <- sum(id == 1)
  N <- length(id)
  
  data_list <- list(
    N = N,
    J = J,
    id = id,
    x = x,
    y = y
  )
  
  fit <- sampling(
    model,
    data = data_list,
    chains = 4,
    iter = 1000,
    warmup = 500,
    seed = 123
  )
  
  post <- rstan::extract(fit)
  
  pred_mu_alpha_mean <- mean(post$mu_alpha)
  rmse_alpha <- sqrt(mean((pred_mu_alpha_mean - true_params$mu_alpha)^2))
  
  pred_mu_beta1_mean <- mean(post$mu_beta1)
  rmse_beta1 <- sqrt(mean((pred_mu_beta1_mean - true_params$mu_beta1)^2))
  
  pred_mu_beta2_mean <- mean(post$mu_beta2)
  rmse_beta2 <- sqrt(mean((pred_mu_beta2_mean - true_params$mu_beta2)^2))
  
  pred_mu_psi_mean <- mean(post$mu_psi)
  rmse_psi <- sqrt(mean((pred_mu_psi_mean - true_params$mu_psi)^2))
  
  rmse_by_size[i,] <- c(rmse_alpha, rmse_beta1, rmse_beta2, rmse_psi, n)
  
  alpha_ci <- quantile(post$mu_alpha, c(0.025, 0.975))
  beta1_ci <- quantile(post$mu_beta1, c(0.025, 0.975))
  beta2_ci <- quantile(post$mu_beta2, c(0.025, 0.975))
  psi_ci <- quantile(post$mu_psi, c(0.025, 0.975))
  
  cover_alpha <- true_params$mu_alpha >= alpha_ci[1] & true_params$mu_alpha <= alpha_ci[2]
  cover_beta1 <- true_params$mu_beta1 >= beta1_ci[1] & true_params$mu_beta1 <= beta1_ci[2]
  cover_beta2 <- true_params$mu_beta2 >= beta2_ci[1] & true_params$mu_beta2 <= beta2_ci[2]
  cover_psi <- true_params$mu_psi >= psi_ci[1] & true_params$mu_psi <= psi_ci[2]
  
  coverage_by_size[i,] <- c(cover_alpha, cover_beta1, cover_beta2, cover_psi, n)
}

rmse_by_size <- as.data.frame(rmse_by_size)
coverage_by_size <- as.data.frame(coverage_by_size)

colnames(rmse_by_size) <- c("alpha", "beta1", "beta2", "psi", "n")
colnames(coverage_by_size) <- c("alpha", "beta1", "beta2", "psi", "n")

plotting_data <- rmse_by_size %>%
  pivot_longer(cols = 1:4, names_to = "param", values_to = "rmse")

params_converge <- plotting_data %>%
  ggplot() +
  geom_line(aes(x = n, y = rmse, color = param)) +
  labs(x = "sample size", color = "parameter")

ggsave("figures/consistency_params.png", plot = params_converge)

coverage_by_size[,1:4] <- lapply(coverage_by_size[,1:4], as.numeric)

coverage_plot_data <- coverage_by_size %>%
  pivot_longer(cols = 1:4, names_to = "param", values_to = "coverage")

coverage_plot <- ggplot(coverage_plot_data) +
  geom_line(aes(x = n, y = coverage, color = param)) +
  labs(x = "sample size", y = "95% CI coverage", color = "parameter") +
  ylim(0,1)

ggsave("figures/consistency_coverage.png", plot = coverage_plot)
