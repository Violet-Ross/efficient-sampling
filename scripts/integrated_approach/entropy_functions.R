compute_entropy <- function(x_star, infection_num) {
  alpha <- post$alpha[, infection_num]
  beta1 <- post$beta1[, infection_num]
  beta2 <- post$beta2[, infection_num]
  psi   <- post$psi[,   infection_num]
  
  mu <- alpha + beta1 * x_star + beta2 * pmax(x_star - psi, 0)
  
  dens <- density(mu)
  pk   <- dens$y / sum(dens$y)
  -sum(pk * log(pk + 1e-10))
}

entropy_maximizer <- function(infection_num, N = 1){
  infx <- all_trajectories %>% filter(index_values == infection_num)
  x_grid <- seq.int(min(infx$time_values), max(infx$time_values), length.out = max(infx$time_values) - min(infx$time_values) + 1)
  
  entropy_vals <- c()
  for(x_val in x_grid){
    entropy_vals <- append(entropy_vals, compute_entropy(x_val, infection_num))
  }
  
  top_indices <- order(entropy_vals, decreasing = TRUE)[1:N]
  x_max_vals <- x_grid[top_indices]
  
  return(x_max_vals)
}

test_points <- function(person_id, x_vals){
  person <- matrix[person_id,]
  infection_start <- all_trajectories[all_trajectories$person_values == person_id, 5][1]
  y_vals <- unlist(matrix[2, (infection_start + x_vals - 1)])
  
  return(y_vals)
}
