# ============================================================================
# utils.R — rank-format conversions, calculation & Inverse-Gamma helpers
# Helpers shared by the samplers, data generators and metrics.
# Part of the BayeSPeer function library; logic unchanged from the original implementation.
# ============================================================================

## --- Rank-format conversions & self/peer diagnostics ----------------------
rank_diff_self_vs_peers <- function(rank_mat) {
  stopifnot(is.matrix(rank_mat),
            nrow(rank_mat) == ncol(rank_mat))
  
  # 1. Copy and blank out the diagonal
  peer_mat <- rank_mat
  diag(peer_mat) <- NA      # leave only peer ranks in each row
  
  # 2. Average peer ranks for each item (NA-aware)
  peer_means <- rowMeans(peer_mat, na.rm = TRUE)
  peer_ranks <- rank(peer_means)
  
  # 3. Difference: self rank – average peer rank
  diff1 <- peer_means - diag(rank_mat) 
  diff2 <- peer_ranks - diag(rank_mat) 
  
  # keep row names (items) if they exist
  if (!is.null(rownames(rank_mat))) names(diff) <- rownames(rank_mat)
  
  return(list(meandiff=diff1, rankdiff=diff2, peer_means=peer_means, peer_ranks=peer_ranks))
}
#
convert_to_rank_lists <- function(rank_matrix) {
  # 1) For each column, sort the ranks and extract the names 
  #    (i.e., the row names of `rank_matrix`).
  #    This returns a character matrix.
  sorted_mat <- apply(rank_matrix, 2, function(col) names(sort(col)))
  
  # 2) Convert each column of that matrix into an element of a list,
  #    so we get a list of length = number of columns.
  lapply(seq_len(ncol(sorted_mat)), function(j) sorted_mat[, j])
}

# # Generate rank lists
# rank_lists <- convert_to_rank_lists(mydat)

convert_to_ranking <- function(ordinal_list) {
  ranking <- seq_along(ordinal_list)  # Assign ranks based on order
  names(ranking) <- ordinal_list      # Assign item names for clarity
  return(ranking[order(names(ranking))]) # Reorder to match original item order
}
# 
# # Example usage
# ordinal_list <- c("Item1", "Item8", "Item3", "Item5", "Item2", "Item6", "Item7", "Item4")
# convert_to_ranking(ordinal_list)

qtruncgamma=function(a,b,p,shape,rate=1)
{qgamma(pgamma(a,shape=shape,rate=rate)+p*(pgamma(b,shape=shape,rate=rate)-pgamma(a,shape=shape,rate=rate)),shape=shape,rate=rate)}
rtruncgamma=function(n,a,b,shape,rate=1)
{u = runif(n, min = 0, max = 1)
x = qtruncgamma(a=a,b=b,p=u,shape=shape,rate=rate)
return(x)}

## --- Calculation helpers --------------------------------------------------
sort_matrix_columns <- function(n, r=mydat$rank) {
  # Initialize W with random normal values
  initW <- W <- matrix(rnorm(n^2, mean = 0, sd = 1), nrow = n, ncol = n)
  
  # Sorting each column of W and storing in initW
  for (i in 1:n) {
    initW[, i] <- sort(W[, i], decreasing = TRUE)[r[, i]]
  }
  
  return(initW)
}


find_sigmas <- function(rhos,n_teams, perc_bias =0.5, randomRho=FALSE){
  
  if (randomRho){
    rhos = runif(n, 0.15, 0.95) 
  }
  
  sig2_beta =sig2_eps = (1/rhos^2 - 1)
  
  return(list(sig2_beta=sig2_beta,
              sig2_eps=sig2_eps))
}


Est_distBeta <- function(mat_rank=mydat$rank, n=n_teams
                         
                         ){
  
  (epsilon=(0.001)/n)  # Small offset
  
  # Generate the equally spaced sequence from epsilon to (1 - epsilon)
  ( mapped_values <- seq(epsilon, 1 - epsilon, length.out = n_teams) )
  
  interval <- (1 - 2*epsilon)/(n-1)
  (sample_values <- rnorm(n, mean=mapped_values, sd=(interval/4)))
  (sample_values <- ifelse(sample_values>1, 1-epsilon, sample_values))
  (sample_values <- ifelse(sample_values<0, epsilon, sample_values))
  
  ( mapped_matrix <- matrix(sample_values[mat_rank], nrow = n, byrow =F) )
  create.omega <- -qnorm(mapped_matrix)
  
  #range(create.omega)
  #range(mydat$weight)
  # (mu.hats <- apply(create.omega, 1, function(x) mean(x, na.rm = TRUE)))
  omega_jj <- diag(create.omega)
  
  omega.selfrm <- create.omega
  diag(omega.selfrm) <- NA
  mu.hats.selfrm <- apply(omega.selfrm, 1, function(x) mean(x, na.rm = TRUE))
  
  # (mubeta.hat <- sum( diag(create.omega)-create.omega )/( n_teams*(n_teams-1) ))
  # (beta.hats <- (mu.hats - mu.hats.selfrm) * n_teams )
  
  beta.hats <- omega_jj - mu.hats.selfrm
  
  return(list(Beta.bar=mean(beta.hats),
              Beta.sd= sd(beta.hats),
              Epsilon.sd= sd(as.numeric(omega.selfrm), na.rm=TRUE))  )
}





## --- Inverse-Gamma hyper-parameter helpers --------------------------------
# Function to estimate a and b parameters of an Inverse Gamma distribution
estimate_inverse_gamma_params <- function(mean_val, perc_1, perc_99) {
  
  # Objective function to minimize
  objective_fn <- function(param) {
    a <- param[1]
    b <- (a - 1) * mean_val
    
    # Compute the 1st and 99th percentiles from the inverse gamma distribution
    perc_1_est <- qinvgamma(0.01, alpha = a, beta = b)
    perc_99_est <- qinvgamma(0.99, alpha = a, beta = b)
    
    # Compute the error between estimated and given percentiles
    error <- (perc_1_est - perc_1)^2 + (perc_99_est - perc_99)^2
    print(param)
    print(error)
    return(error)
  }
  
  # Initial guess for a
  init_a <- 2.1
  
  # Optimization
  result <- optim(init_a, objective_fn, method = "L-BFGS-B", lower = 1.01)
  a_est <- result$par
  b_est <- (a_est - 1) * mean_val
  
  # Return estimated parameters
  return(list(a = a_est, b = b_est))
}


rinvgamma_trunc <- function(n, shape, scale, lower, upper) {
  # Compute the CDF at the truncation points
  F_lower <- pinvgamma(lower, alpha = shape, beta = scale)
  F_upper <- pinvgamma(upper, alpha = shape, beta = scale)
  
  # Generate uniform random variables between F_lower and F_upper
  U <- runif(n, min = F_lower, max = F_upper)
  
  # Apply the inverse CDF (quantile function) to get samples
  samples <- qinvgamma(U, alpha = shape, beta = scale)
  
  return(samples)
}

find_shape_rate <- function(mean, sd){ #for Inverse-Gamma
  
  var <- sd^2
  alpha <- (mean^2) / var + 2
  
  # Calculate beta using the formula for mean
  beta <- mean * (alpha - 1)
  
  return(list(shape = alpha, rate = beta))
  
}

#plot( y=(1/pgamma(seq(0.001, 10, 0.001),shape=0.01, rate=0.01 )),  x=seq(0.001, 10, 0.001))
# find_shape_rate(mean=0.2, sd=0.05)
# find_shape_rate(mean=2, sd=1.5)
#
# 
# quantile(1/rgamma(20000,
#                   shape=18,
#                   rate=3.4 ), probs = c(0, 0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975))
# quantile(1/rgamma(20000,
#                   shape=3.78,
#                   rate=5.56 ), probs = c(0, 0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975))
# 
# mean(1/rgamma(20000,
#               shape=8.25,
#               rate=1.8125 ) )
# sd(1/rgamma(20000,
#               shape=8.25,
#               rate=1.8125 ) )
# 
# mean(1/rgamma(20000,
#               shape=4.25,
#               rate=9.75 ) )
# sd(1/rgamma(20000,
#              shape=4.25,
#              rate=9.75 ) )

