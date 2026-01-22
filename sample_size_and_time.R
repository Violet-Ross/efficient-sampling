library(tidyverse)
library(segmented)

## Function to Generate Trajectories
generate_traj <- function(num_traj){
  t <- c() # vector of times
  vol <- c() # vector of viral load
  trajec_num <- c() # vector of trajectory indices
  
  for(i in 1:num_traj){ 
    # generate t_0 and append to list
    t_0 <- runif(1, 2.5, 3.5)
    t <- append(t, t_0)
    
    candidate_w_1 <- Inf
    while(candidate_w_1 > 2.5){
      candidate_w_1 <- rgamma(1, 1.5, 1)
    }
    t <- append(t, t_0 + candidate_w_1 + 0.5)
    w_2 <- runif(1, 4, 9)
    t <- append(t, t_0 + candidate_w_1 + 0.5 + w_2)
    
    v_peak <- runif(1, 7, 11)
    vol <- append(vol, c(3, v_peak, 6))
    
    trajec_num <- append(trajec_num, c(i, i, i))
  } 
  
  df <- as.data.frame(cbind(t, vol, trajec_num))
  
  return(df)
}

## Function to return actual y values for x values
get_y <- function(traj, x){
  m_1 <-(traj$vol[[2]] - traj$vol[[1]]) / (traj$t[[2]] - traj$t[[1]]) # slope of first line seg
  b_1 <- traj$vol[[2]] - (m_1 * traj$t[[2]]) # y-int of first line seg
  
  m_2 <-(traj$vol[[3]] - traj$vol[[2]]) / (traj$t[[3]] - traj$t[[2]]) # slope of second line seg
  b_2 <- traj$vol[[3]] - (m_2 * traj$t[[3]]) # y-int of second line seg
  
  y <- c()
  for(i in 1:length(x)){
    if(x[i] <= traj$t[[2]]){
      y_i <- m_1*x[i] + b_1
    }
    else{
      y_i <- m_2*x[i] + b_2
    }
    y <- c(y, y_i)
  }
  
  return(y)
}

## Function to sample from trajectory
sample_traj <- function(traj, num_samples){
  x <- runif(num_samples, min = traj$t[[1]], max = traj$t[[3]]) # generate u.a.r. x coords
  
  y <- get_y(traj, x)
  
  sampled_pts <- data.frame(x = x, y = y)
  
  return(sampled_pts)
}

## Simulation
max_n = 20 
num_iter = 100

all_mse <- c()

for(n in 2:max_n){
  print(n)
  sink(file = "/dev/null", append = FALSE, type = "output")
  mse_for_n <- c()
  for(i in 1:num_iter){    
    traj <- generate_traj(1)
    sampled_pts <- sample_traj(traj, n)
    fit.lm <- lm(data = sampled_pts, y ~ x)
    fit.seg  <- suppressWarnings(selgmented(fit.lm, Kmax = 1, type = "bic", check.dslope=FALSE))
    
    to_predict <- seq(3, max(traj$t), by = 0.25)
    fit.data <- data.frame(x = to_predict, y = predict(fit.seg, data.frame(x = to_predict)))
    actual <- get_y(traj, to_predict)
    
    mse <- (1/length(to_predict)) * sum((actual - fit.data$y)^2)
    mse_for_n <- c(mse_for_n, mse)
  }
  all_mse <- c(all_mse, mean(mse_for_n))
  sink()
}

plot1 <- ggplot() +
  geom_point(aes(x = 2:n, y = all_mse), size = 5) +
  labs(y = "mse", x = "number of samples") +
  theme_minimal() +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 25))
  

ggsave("/Users/violetross/Desktop/FTAS/Sampling/efficient-sampling/efficient-sampling/figures/sample_size.png",
       plot = plot1, width = 20, height = 9)

## Why are my results so weird? AKA modified simulation
max_n = 20 
num_iter = 20
should_break <- FALSE

all_mse <- c()

for(n in 2:max_n){
  print(n)
  sink(file = "/dev/null", append = FALSE, type = "output")
  mse_for_n <- c()
  for(i in 1:num_iter){    traj <- generate_traj(1)
  sampled_pts <- sample_traj(traj, n)
  fit.lm <- lm(data = sampled_pts, y ~ x)
  fit.seg  <- suppressWarnings(selgmented(fit.lm, Kmax = 1, type = "bic", check.dslope=FALSE))
  
  to_predict <- seq(3, max(traj$t), by = 0.25)
  fit.data <- data.frame(x = to_predict, y = predict(fit.seg, data.frame(x = to_predict)))
  actual <- get_y(traj, to_predict)
  
  mse <- (1/length(to_predict)) * sum((actual - fit.data$y)^2)
  if(mse > 50){
    plot <- ggplot() +
      geom_point(aes(x = to_predict, y = fit.data$y), size = 5) +
      geom_line(aes(x = to_predict, y = actual), linewidth = 3)
    sink()
    should_break <- TRUE
    break
  }
  mse_for_n <- c(mse_for_n, mse)
  }
  if(should_break == TRUE){
    break
  }
  all_mse <- c(all_mse, mean(mse_for_n))
  sink()
}

plot2 <- plot +
  labs(y = "viral RNA", x = "t") +
  theme_minimal() +
  theme(axis.text = element_text(size = 14), axis.title = element_text(size = 25))

ggsave("/Users/violetross/Desktop/FTAS/Sampling/efficient-sampling/efficient-sampling/figures/sample_time.png",
       plot = plot2, width = 20, height = 9)
