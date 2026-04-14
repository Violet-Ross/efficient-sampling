library(tidyverse)

setwd("/Users/violetross/Desktop/FTAS/efficient-sampling")

M_0 <- 10 # number of individuals to use for the first round of training
# M_0 is NOT the same as number of infections, since many of these individuals will
# not have an infection
M <- 4 # number of new infections to fit to at each iteration
N_add <- 2 # number of additional samples to take from each infection

matrix <- read_csv("throughput/synthetic_matrix.csv")
true_params <- read_csv("throughput/integrated_true_population_params.csv")

source("scripts/integrated_approach/initial_infx.R")
source("scripts/integrated_approach/find_infx.R")

to_bind <- partial_trajectories %>% 
  mutate(index_values = index_values + max(init_trajectories$index_values)) %>%
  mutate(person_values = person_values + max(init_trajectories$person_values))

all_trajectories <- rbind(init_trajectories, to_bind)

source("scripts/integrated_approach/viral_bhm.R")

post  <- extract(fit)
source("scripts/integrated_approach/entropy_functions.R")

i <- max(init_trajectories$index_values) + 1
rmse_by_iteration <- data.frame(rmse_alpha = numeric(),
                                rmse_beta1 = numeric(),
                                rmse_beta2 = numeric(),
                                rmse_psi = numeric())

fitting_data <- all_trajectories
while(i < max(all_trajectories$index_values)){
  next_infx <- all_trajectories %>% filter(index_values %in% i : (i + M - 1))
  for(infx_id in unique(next_infx$index_values)){
    print(infx_id)
    samples_x <- entropy_maximizer(infx_id, N_add)
    person_id <- all_trajectories %>% 
      filter(index_values == infx_id) %>% 
      pull(4) %>%
      first()
    samples_y <- test_points(person_id, samples_x)
    for(j in 1:length(samples_x)){
      new_row <- data.frame(index_values = infx_id, 
                            time_values = samples_x[j], 
                            traj_values = samples_y[j], 
                            person_values = NA,
                            date_values = NA)
      fitting_data <- rbind(fitting_data, new_row)
    }
  }
  fitting_data <- arrange(fitting_data, index_values, time_values)
  source("scripts/integrated_approach/viral_bhm.R")
  
  post  <- extract(fit)
  alpha_rmse <- sqrt(mean((post$mu_alpha - true_params$mu_alpha)^2))
  beta1_rmse <- sqrt(mean((post$mu_beta1 - true_params$mu_beta1)^2))
  beta2_rmse <- sqrt(mean((post$mu_beta2 - true_params$mu_beta2)^2))
  psi_rmse   <- sqrt(mean((post$mu_psi   - true_params$mu_psi)^2))
  this_rmse <- c(alpha_rmse, beta1_rmse, beta2_rmse, psi_rmse)

  rmse_by_iteration <- rbind(rmse_by_iteration, this_rmse)
  
  i = i + M
}

colnames(rmse_by_iteration) <- c("alpha", "beta1", "beta2", "psi")

pivoted_rmse <- rmse_by_iteration %>%
  pivot_longer(cols = 1:4, names_to = "variable", values_to = "rmse") %>%
  mutate(iter = rep(1:nrow(rmse_by_iteration), each = 4))

pivoted_rmse %>%
  ggplot() +
    geom_line(aes(x = iter, y = rmse, color = variable)) +
    labs(x = "iteration", y = "rmse", title = "Experiment 1:\nM_0 = 10, M = 4, N_add = 2")

ggsave("figures/rmse_by_iter_1.png")
