library(knitr)
library(kableExtra)

# Define function that computes the KL divergence of two distributions
# from one vector of samples from each distribution
kl_divergence <- function(true_dist_samples, posterior_samples){
  # Kernel density estimates
  dx <- density(true_dist_samples, n = 2048)
  dy <- density(posterior_samples, n = 2048)
  
  # Interpolate Q density onto P grid
  q_interp <- approx(dy$x, dy$y, xout = dx$x, rule = 2)$y
  
  # Avoid zeros (important!)
  eps <- 1e-12
  p <- dx$y + eps
  q <- q_interp + eps
  
  # Normalize (ensures they integrate to 1 numerically)
  p <- p / sum(p)
  q <- q / sum(q)
  
  # KL divergence
  kl_pq <- sum(p * log(p / q))
  kl_pq
}

# read in the fit from viral_bhm.R
fit <- readRDS("stan_fit.rds")
post <- rstan::extract(fit)

# read in the true params used to generate the synthetic data
true_params <- read_csv("true_population_params.csv")

# get actual distributions
actual_alpha <- rnorm(1000000, true_params$mu_alpha, true_params$sigma_alpha)
actual_beta1 <- rnorm(1000000, true_params$mu_beta1, true_params$sigma_beta1)
actual_beta2 <- rnorm(1000000, true_params$mu_beta2, true_params$sigma_beta2)
actual_psi <- rnorm(1000000, true_params$mu_psi, true_params$sigma_psi)

# compute divergence for each parameter's distribution
true_dists <- rbind(actual_alpha, actual_beta1, actual_beta2, actual_psi)
posterior_dists <- rbind(c(post$alpha), c(post$beta1), c(post$beta2), c(post$psi))
divergences <- c()
for(i in 1:4){
  divergences <- c(divergences, kl_divergence(true_dists[i,], posterior_dists[i,]))
}

# put the divergence information into a table and save the table as a .png
divergences_by_param <- data.frame(parameter = c("alpha", "beta1", "beta2", "psi"), divergences)

divergence_table <- kable(divergences_by_param, caption = "KL divergence of true and posterior distributions") %>%
  kable_classic(full_width = F) %>%
  column_spec(1:2, width = "150px")

save_kable(divergence_table, 
           zoom = 4,
           density = 300,
           "divergence_table.png")
  
