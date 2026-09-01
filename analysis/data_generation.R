# ============================================================================
# data_generation.R — synthetic data generators
# gen_data_dir / gen_data_mix / gen_data_t / gen_data_gamma / gen_data_tri.
# Part of the BayeSPeer function library.
# ============================================================================

gen_data_dir <- function(
    n = 5,   # number of teams
    mu_beta,
    sigma_beta,
    sigma_epsilon
) {
  regenerate <- TRUE  # Flag to control regeneration loop
  result <- NULL
  
  while (regenerate) {
    # Generate data
    mu <- rnorm(n, mean = 0, sd = 1)  # Global team effects
    true_rank <- rank(-mu)
    
    beta <- rnorm(n, mean = mu_beta, sd = sigma_beta)
    epsilon <- matrix(rnorm(n * n, mean = 0, sd = sigma_epsilon), nrow = n, ncol = n)
    
    diag_effects <- diag(beta)
    omega <- outer(mu, rep(1, n)) + diag_effects + epsilon  # Fully vectorized
    
    
    omega_rank <- apply(-omega, 2, rank)
    omega_var <- matrixStats::colVars(omega) 
    rhosq <- 1 / omega_var
    
    rownames(omega_rank) <- paste0("Item", 1:n)
    colnames(omega_rank) <- paste0("Ranker", 1:n)
    
    # Check for rows with identical rank values
    identical_rows <- (nrow(unique(omega_rank))==1)
    
    if (!identical_rows) {
      regenerate <- FALSE  # Exit the loop if no identical rows
      result <- list(
        weight = omega,
        rank = omega_rank,
        rank_list = convert_to_rank_lists(omega_rank),
        true_mu = mu,
        true_rank = true_rank,
        var_w = omega_var,
        rhosq = rhosq,
        true_beta = beta
      )
    }
  }
  
  return(result)
}


gen_data_mix <- function(
    n = 5,   # number of teams
    mu_beta,
    sigma_beta
    # sigma_epsilon removed — now generated internally per ranker
) {
  regenerate <- TRUE  # Flag to control regeneration loop
  result <- NULL
  
  while (regenerate) {
    # Generate data
    mu <- rnorm(n, mean = 0, sd = 1)  # Global team effects
    true_rank <- rank(-mu)
    
    beta <- rnorm(n, mean = mu_beta, sd = sigma_beta)
    
    # Each ranker gets their own rho drawn from U(0.15, 0.95)
    rhos <- runif(n, min = 0.15, max = 0.95)
    sigma2_epsilons <- (1/rhos^2 - 1)   # variances: sigma_j^2 = 1/rho_j^2 - 1
    
    # Build epsilon matrix column-by-column, one sigma per ranker (column)
    epsilon <- sapply(sigma2_epsilons, function(s2) rnorm(n, mean = 0, sd = sqrt(s2)))
    
    diag_effects <- diag(beta)
    omega <- outer(mu, rep(1, n)) + diag_effects + epsilon  # Fully vectorized
    
    omega_rank <- apply(-omega, 2, rank)
    omega_var <- matrixStats::colVars(omega) 
   # rhosq <- 1 / omega_var #maybe wrong not using the theoretical formula and esp with n=5 very noisy
    
    rownames(omega_rank) <- paste0("Item", 1:n)
    colnames(omega_rank) <- paste0("Ranker", 1:n)
    
    # Check for rows with identical rank values
    identical_rows <- (nrow(unique(omega_rank)) == 1)
    
    if (!identical_rows) {
      regenerate <- FALSE  # Exit the loop if no identical rows
      result <- list(
        weight = omega,
        rank = omega_rank,
        rank_list = convert_to_rank_lists(omega_rank),
        true_mu = mu,
        true_rank = true_rank,
        true_rhos=rhos,
        var_w = omega_var,
       # rhosq = rhosq,
        true_beta = beta,
        sigma2_epsilons = sigma2_epsilons  # ranker noise variances, for inspection
      )
    }
  }
  
  return(result)
}
#------------------------
# Robustness
gen_data_t <- function(
    n = 5,   # number of teams
    mu_beta,
    sigma_beta,
    sigma_epsilon,
    df
) {
  regenerate <- TRUE  # Flag to control regeneration loop
  result <- NULL
  
  while (regenerate) {
    # Generate data
    mu <- rt(n, df = df) * sqrt((df - 2) / df) # Global team effects, rescaled to unit variance
    true_rank <- rank(-mu)
    
    beta <- rnorm(n, mean = mu_beta, sd = sigma_beta)
    epsilon <- matrix(rnorm(n * n, mean = 0, sd = sigma_epsilon), nrow = n, ncol = n)
    
    diag_effects <- diag(beta)
    omega <- outer(mu, rep(1, n)) + diag_effects + epsilon  # Fully vectorized
    
    omega_rank <- apply(-omega, 2, rank)
    omega_var <- apply(omega, 2, var)
    rhosq <- 1 / omega_var
    
    rownames(omega_rank) <- paste0("Item", 1:n)
    colnames(omega_rank) <- paste0("Ranker", 1:n)
    
    # Check for rows with identical rank values
    identical_rows <- (nrow(unique(omega_rank))==1)
    
    if (!identical_rows) {
      regenerate <- FALSE  # Exit the loop if no identical rows
      result <- list(
        weight = omega,
        rank = omega_rank,
        rank_list = convert_to_rank_lists(omega_rank),
        true_mu = mu,
        true_rank = true_rank,
        var_w = omega_var,
        rhosq = rhosq,
        true_beta = beta
      )
    }
  }
  
  return(result)
}


gen_data_gamma <- function(
    n = 5,   # number of teams
    mu_beta,
    sigma_beta,
    sigma_epsilon,
    kappa
) {
  regenerate <- TRUE  # Flag to control regeneration loop
  result <- NULL
  
  while (regenerate) {
    # Generate data
    mu <- rgamma(n, shape = kappa, rate = sqrt(kappa))# Global team effects
    true_rank <- rank(-mu)
    
    beta <- rnorm(n, mean = mu_beta, sd = sigma_beta)
    epsilon <- matrix(rnorm(n * n, mean = 0, sd = sigma_epsilon), nrow = n, ncol = n)
    
    diag_effects <- diag(beta)
    omega <- outer(mu, rep(1, n)) + diag_effects + epsilon  # Fully vectorized
    
    omega_rank <- apply(-omega, 2, rank)
    omega_var <- apply(omega, 2, var)
    rhosq <- 1 / omega_var
    
    rownames(omega_rank) <- paste0("Item", 1:n)
    colnames(omega_rank) <- paste0("Ranker", 1:n)
    
    # Check for rows with identical rank values
    identical_rows <- (nrow(unique(omega_rank))==1)
    
    if (!identical_rows) {
      regenerate <- FALSE  # Exit the loop if no identical rows
      result <- list(
        weight = omega,
        rank = omega_rank,
        rank_list = convert_to_rank_lists(omega_rank),
        true_mu = mu,
        true_rank = true_rank,
        var_w = omega_var,
        rhosq = rhosq,
        true_beta = beta
      )
    }
  }
  
  return(result)
}


# Archieved
# gen_data_dir <- function(# Parameters
#   n= 10,   # number of teams
#   # Hyper-parameters
#   mu_beta,
#   # beta = rep(c(-1, 1), each=5),
#   sigma_beta,
#   sigma_epsilon#,
#   # perc_unbiased=0
# ){
#   
#   mu <- rnorm(n, mean = 0, sd = 1)                  # global team effects
#   true_rank <- rank(-mu)
#   # model components
#   beta <- rnorm(n, mean = mu_beta, sd = sigma_beta)
#   
#   # not everybody has bias
#   # n_unbiased <- round(n*perc_unbiased)
#   # beta[sample(n_unbiased)] <- 0
#   
#   epsilon <- matrix(rnorm(n * n, mean = 0, sd = sigma_epsilon), nrow = n, ncol = n)
#   
#   # Initialize omega matrix
#   omega <- rank <- matrix(NA, nrow = n, ncol = n)
#   
#   # Compute omega for each observation
#   for (i in 1:n) {
#     for (j in 1:n) {
#       
#       omega[i, j] <- mu[i] + (i == j) * beta[j] + epsilon[i,j]
#     }
#   }
#   omega_rank <- apply(-omega, 2, rank)
#   omega_var <- apply(omega, 2, var)
#   rhosq <- 1/omega_var
#   
#   return(list(weight=omega,
#               rank=omega_rank, 
#               true_mu = mu, 
#               true_rank=true_rank,
#               var_w = omega_var, 
#               rhosq=rhosq, 
#               true_beta=beta))
# }

gen_data_tri <- function(
    n = 5,   # number of teams
    mu_beta_vec = c(over = 1, unbiased = 0, under = -1), # Default example offsets
    sigma_beta,
    sigma_epsilon
) {
  regenerate <- TRUE  
  result <- NULL
  
  while (regenerate) {
    # 1. Global team effects
    mu <- rnorm(n, mean = 0, sd = 1)  
    true_rank <- rank(-mu)
    
    # 2. Assign Ranker Types (Tri-modal)
    # 50% Over-confident (index 1), 25% Unbiased (index 2), 25% Under-confident (index 3)
    types <- sample(1:3, size = n, replace = TRUE, prob = c(0.50, 0.25, 0.25))
    
    # 3. Generate beta based on assigned type
    # We map each ranker's type to the corresponding mean in mu_beta_vec
    beta_means <- mu_beta_vec[types]
    beta <- rnorm(n, mean = beta_means, sd = sigma_beta)
    
    # 4. Generate noise and combine
    epsilon <- matrix(rnorm(n * n, mean = 0, sd = sigma_epsilon), nrow = n, ncol = n)
    diag_effects <- diag(beta)
    omega <- outer(mu, rep(1, n)) + diag_effects + epsilon  
    
    # 5. Ranking and Stats
    omega_rank <- apply(-omega, 2, rank)
    omega_var <- matrixStats::colVars(omega) 
    rhosq <- 1 / omega_var
    
    rownames(omega_rank) <- paste0("Item", 1:n)
    colnames(omega_rank) <- paste0("Ranker", 1:n)
    
    # Check for validity (identical rows check)
    identical_rows <- (nrow(unique(omega_rank)) == 1)
    
    if (!identical_rows) {
      regenerate <- FALSE 
      result <- list(
        weight = omega,
        rank = omega_rank,
        rank_list = convert_to_rank_lists(omega_rank),
        true_mu = mu,
        true_rank = true_rank,
        var_w = omega_var,
        rhosq = rhosq,
        true_beta = beta,
        ranker_types = types # Added this so you can track who was what
      )
    }
  }
  
  return(result)
}
