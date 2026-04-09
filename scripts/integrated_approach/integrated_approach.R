library(tidyverse)

setwd("/Users/violetross/Desktop/FTAS/efficient-sampling")

M_0 <- 100 # number of individuals to use for the first round of training
# M_0 is NOT the same as number of infections, since many of these individuals will
# not have an infection
M <- 4 # number of new infections to fit to at each iteration
N <- 2 # number of additional samples to take from each infection

matrix <- read_csv("throughput/synthetic_matrix.csv")

source("scripts/integrated_approach/initial_infx.R")
source("scripts/integrated_approach/find_infx.R")

to_bind <- pooled_trajectories %>% 
  mutate(index_values = index_values + max(init_trajectories$index_values)) %>%
  mutate(person_values = person_values + max(init_trajectories$person_values))

all_trajectories <- rbind(init_trajectories, to_bind)

source("scripts/integrated_approach/viral_bhm.R")

post  <- extract(fit)
source("scripts/integrated_approach/entropy_functions.R")

i <- max(init_trajectories$index_values) + 1
rmse_list <- c()

fitting_data <- all_trajectories
while(i < 38){
  next_infx <- all_trajectories %>% filter(index_values %in% i : (i + M - 1))
  for(infx_id in unique(next_infx$index_values)){
    samples_x <- entropy_maximizer(infx_id, N)
    person_id <- all_trajectories[all_trajectories$index_values == infx_id, 4][1]
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
  rmse_list <- c(rmse_list, this_rmse)
}
