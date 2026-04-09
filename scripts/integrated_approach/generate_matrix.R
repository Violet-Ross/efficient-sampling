num_people <- 500
num_days <- 500

prev <- 0.3
sigma_epsilon = 0.1 # amount of random noise on observations
min_x <- 0
max_x <- 50

# initialize matrix of 0s
infection_matrix <- matrix(0, nrow = num_people, ncol = num_days)

# determine which rows will be infected
is_infx <- rbinom(num_people, size = 1, prev)

# population level parameters (drawn from priors)
mu_alpha <- rnorm(n = 1, mean = 2, sd = 0.5)
mu_beta1 <- rnorm(1, 2, 0.1)
mu_beta2 <- rnorm(1, -4, 0.1)
mu_psi <- rnorm(1, 5, 0.5)

sigma_alpha <- abs(rnorm(1, 0, sd = 0.1))
sigma_beta1 <- abs(rnorm(1, 0, sd = 0.1))
sigma_beta2 <- abs(rnorm(1, 0, sd = 0.1))
sigma_psi <- abs(rnorm(1, 0, sd = 0.1))

for(row_ind in 1:nrow(infection_matrix)){
  if(is_infx[row_ind] == 1){
    # individual level parameters
    alpha <- rnorm(1, mu_alpha, sigma_alpha)
    beta1 <- rnorm(1, mu_beta1, sigma_beta1)
    beta2 <- rnorm(1, mu_beta2, sigma_beta2)
    psi <- rnorm(1, mu_psi, sigma_psi)
    
    x_vals = seq(from = min_x, to = max_x, length.out = max_x - min_x)
    
    # observation level
    y_vals <- c()
    for(x in x_vals){
      obs_mu = alpha + (beta1 * x) + (beta2 * (x - psi) * (x > psi))
      y <- rnorm(1, obs_mu, sigma_epsilon)
      y_vals <- c(y_vals, y)
    } 
    pos_y_vals <- y_vals[y_vals > 0]
    infx_length <- length(pos_y_vals)
    
    start_day <- sample(1:(num_days - infx_length), 1)
    
    infection_matrix[row_ind, start_day:(start_day + infx_length - 1)] <- pos_y_vals
  }
}

write_csv(data.frame(infection_matrix), "throughput/synthetic_matrix.csv")
