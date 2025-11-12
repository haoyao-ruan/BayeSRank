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


# Robustness
gen_data_t <- function(
    n = 5,   # number of teams
    mu_beta,
    sigma_beta,
    sigma_epsilon
) {
  regenerate <- TRUE  # Flag to control regeneration loop
  result <- NULL
  
  while (regenerate) {
    # Generate data
    df <- 5
    mu <- rt(n, df = df) * sigma_beta / sqrt(df / (df - 2)) + mu_beta # Global team effects
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


gen_data_absnorm <- function(
    n = 5,   # number of teams
    mu_beta,
    sigma_beta,
    sigma_epsilon
) {
  regenerate <- TRUE  # Flag to control regeneration loop
  result <- NULL
  
  while (regenerate) {
    # Generate data
    mu <- abs(rnorm(n, mean = 0, sd = 1) ) # Global team effects
    true_rank <- rank(-mu)
    
    beta <- rnorm(n, mean = mu_beta, sd = sigma_beta)
    epsilon <- matrix(rnorm(n * n, mean = 0, sd = sigma_epsilon), nrow = n, ncol = n)
    
    omega <- rank <- matrix(NA, nrow = n, ncol = n)
    for (i in 1:n) {
      for (j in 1:n) {
        omega[i, j] <- mu[i] + (i == j) * beta[j] + epsilon[i, j]
      }
    }
    
    omega_rank <- apply(-omega, 2, rank)
    omega_var <- apply(omega, 2, var)
    rhosq <- 1 / omega_var
    
    # Check for rows with identical rank values
    identical_rows <- any(apply(omega_rank, 1, function(row) length(unique(row)) == 1))
    
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
