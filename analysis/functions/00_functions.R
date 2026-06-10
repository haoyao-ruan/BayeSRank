#-------------------------------------
library(MCMCpack)
library(extraDistr)

library(TopKLists)
library("expm")  # For matrix calculations used in the geometric mean

library(tidyverse)
library(reshape2)
library(matrixStats) #colVar

#
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
#------------------Gibbs Sampler Application-------------------
BayeSRank.pre <- function(n, r, M, burnin, 
                      # verbose,
                      W,             # W=initial Omega matrix
                      mu0,
                      mubeta0,
                      beta0,
                      Var_beta0,
                      Var_epsilon0, minbeta, maxbeta, 
                      a,b,c,d
){
  
  # Dimensions of the rank matrix
  N <- n  # Number of teams (rows of `r`)
  J <- n  # Number of rankers (columns of `r`)
  
  # Initialize matrices to store Mu, Beta, and Sigma_s2
  Mu <- Beta <- matrix(1, N, M)
  Var_epsilon <- matrix(1, N, M)
  Mu_beta <- Var_beta <- matrix(1, 1, M) # DO NOT USE 0 coz NAs will be produced at 2nd itration
  
  # Initial values
  Mu[, 1] <- mu0
  Mu_beta[1, 1] <- mubeta0
  Var_beta[1, 1] <- Var_beta0
  
  beta0 <- rnorm(J, mubeta0, sd=sqrt(Var_beta0) )
  Beta[, 1] <- beta0
  
  Var_epsilon[, 1] <- Var_epsilon0
  
  for (m in 2:M) {
    
    # if (m %% verbose == 0) { 
    #   print(m)  # Print progress
    # }
    
    #   # 1.Update Mu
    
    for (i in 1:N) { #for each team in ROWS
      
      indBeta <- Beta[, m-1]
      indBeta[-i] <-0 
      
      g <- sum( (W[i, ] - indBeta) / Var_epsilon[, m-1] ) #numerator g
      h <- sum(1 / Var_epsilon[, m-1]) + 1                #denominator h
      
      Mu[i, m] <- rnorm(1, mean = (g/h), sd = sqrt(1/h) ) 
      
    }
    
    Mu[,m] <- Mu[,m]- mean(Mu[,m])
    
    #   # 2.Update Beta
    
    for (j in 1:J) { #for each ranker in COLUMNS
      
      g <- ( (W[j, j] - Mu[j, m]) / Var_epsilon[j, m-1] ) + (Mu_beta[m-1]/Var_beta[m-1]) #numerator g
      h <- (1 / Var_epsilon[j, m-1]) + (1/Var_beta[m-1])              #denominator h
      
      Beta[j, m] <- rnorm(1, (g/h), sqrt(1/h))
      
      # Beta[j, m] <- truncnorm::rtruncnorm(1, a = minbeta*2, b = maxbeta*2,
      #                                     mean = (g/h),
      #                                     sd = sqrt(1/h))
    } 
    
    
    #   # 3.Update mu_beta
    
    # g <- sum(Beta[,m]/Var_beta[m-1])                                      #numerator g
    # h <- (1 / Var_beta[m-1]) + (1/tau)^2 #denominator h
    # 
    # Mu_beta[m] <- rnorm(1, (g/h), sqrt(1/h))
    
    Mu_beta[m] <- truncnorm::rtruncnorm(1, a = minbeta, b = maxbeta,
                                        mean = mean(Beta[, m]),
                                        sd = sqrt(Var_beta[m-1]/J) )
    
    #   # 4.Update Var_beta
    # Var_beta[m] <-  1/rtruncgamma(1, a=0.25, b=Inf,
    #                               shape = a+ (J/2),
    #                               rate = b+ sum( (Beta[, m]-Mu_beta[,m])^2 )/2 )
     Var_beta[m] <-  1/rgamma(1,
                             shape = a+ (J/2),
                             rate = b+ sum( (Beta[, m]-Mu_beta[,m])^2 )/2 )


    #   # 5.Update Var_epsilon
    
    for (j in 1:J) {
      
      indBeta <- Beta[, m]
      indBeta[-j] <-0 
      
      Var_epsilon[j, m] <- 1/rgamma(1,
                                    shape = c+ (J/2),
                                    rate = d+ sum( (W[, j] - Mu[,m] - indBeta)^2 )/2 )

      # Var_epsilon[j, m] <- 1/rtruncgamma(1,a=(1/400), b=Inf,
      #                               shape = c+ (J/2),
      #                               rate = d+ sum( (W[, j] - Mu[,m] - indBeta)^2 )/2 )
       
 
    }
    
    
    #   # 6. Update W
    
    for (j in 1:J) {
      for (i in 1:N) {
        
        if (r[i, j] == 1) {
          upper <- 100
          lower <- W[which(r[, j] == 2), j]
          
        } else if (r[i,j] == n) {
          upper <- W[which(r[, j] == n - 1), j]
          lower <- -100
          
        } else {
          
          lower <- W[which(r[, j] == r[i,j] + 1), j]
          upper <- W[which(r[, j] == r[i,j] - 1), j]
          
        }
        W[i,j] <- truncnorm::rtruncnorm(1, a = lower, b = upper,
                                        mean = Mu[i, m] + Beta[j]*(i==j),
                                        sd = sqrt( Var_epsilon[j, m] ))
      }
    }
  }
  
  drop = 1:burnin
  post.mean.mus <- apply(Mu[,-drop], 1, mean)
  
  post.mean.betas <- apply(Beta[,-drop], 1, mean)
  post.median.Var_epsilon <- apply(Var_epsilon[,-drop], 1, median)
  
  post.mean.mu_beta <- mean(Mu_beta[,-drop])
  post.median.Var_beta <- median(Var_beta[,-drop])
  
  
  out <- list(mu = Mu, 
    Beta_draws = Beta,
    muBeta_draws = Mu_beta,
    varBeta_draws = Var_beta,
    varEpsilon_draws = Var_epsilon,
    #W=W,
    post_mean_mus = post.mean.mus,
    post_mean_betas = post.mean.betas, 
    post_mean_mubeta = post.mean.mu_beta, 
    post_median_Varbeta = post.median.Var_beta,
    post_median_VarEps = post.median.Var_epsilon,
    post_rank = rank(-post.mean.mus)
  )
  
  return(out)
}  

#------------------Gibbs Sampler Simulation only report estimates-------------
BayeSRank <- function(n, r, M, burnin, 
                          # verbose,
                          W,             # W=initial Omega matrix
                          mu0,
                          mubeta0,
                          beta0,
                          Var_beta0,
                          Var_epsilon0, minbeta, maxbeta, 
                          a,b,c,d,
                          sig2_lower, sig2_upper
){
  
  # Dimensions of the rank matrix
  N <- n  # Number of teams (rows of `r`)
  J <- n  # Number of rankers (columns of `r`)
  
  # Initialize matrices to store Mu, Beta, and Sigma_s2
  Mu <- Beta <- matrix(1, N, M)
  Var_epsilon <- matrix(1, N, M)
  Mu_beta <- Var_beta <- matrix(1, 1, M) # DO NOT USE 0 coz NAs will be produced at 2nd itration
  
  # Initial values
  Mu[, 1] <- mu0
  Mu_beta[1, 1] <- mubeta0
  Var_beta[1, 1] <- Var_beta0
  
  beta0 <- rnorm(J, mubeta0, sd=sqrt(Var_beta0) )
  Beta[, 1] <- beta0
  
  Var_epsilon[, 1] <- Var_epsilon0
  
  for (m in 2:M) {
    
    # if (m %% verbose == 0) { 
    #   print(m)  # Print progress
    # }
    
    #   # 1.Update Mu
    
    for (i in 1:N) { #for each team in ROWS
      
      indBeta <- Beta[, m-1]
      indBeta[-i] <-0 
      
      g <- sum( (W[i, ] - indBeta) / Var_epsilon[, m-1] ) #numerator g
      h <- sum(1 / Var_epsilon[, m-1]) + 1                #denominator h
      
      Mu[i, m] <- rnorm(1, mean = (g/h), sd = sqrt(1/h) ) 
      
    }
    
    Mu[,m] <- Mu[,m]- mean(Mu[,m])
    
    #   # 2.Update Beta
    
    for (j in 1:J) { #for each ranker in COLUMNS
      
      g <- ( (W[j, j] - Mu[j, m]) / Var_epsilon[j, m-1] ) + (Mu_beta[m-1]/Var_beta[m-1]) #numerator g
      h <- (1 / Var_epsilon[j, m-1]) + (1/Var_beta[m-1])              #denominator h
      
      Beta[j, m] <- rnorm(1, (g/h), sqrt(1/h))
    } 
    
    
    #   # 3.Update mu_beta
    Mu_beta[m] <- truncnorm::rtruncnorm(1, a = minbeta, b = maxbeta,
                                        mean = mean(Beta[, m]),
                                        sd = sqrt(Var_beta[m-1]/J) )
    
    #   # 4.Update Var_beta
    
    Var_beta[m] <-  1/rgamma(1,
                             shape = a+ (J/2),
                             rate = b+ sum( (Beta[, m]-Mu_beta[,m])^2 )/2 )
    
    
    #   # 5.Update Var_epsilon
    
    for (j in 1:J) {
      
      indBeta <- Beta[, m]
      indBeta[-j] <-0 
      
      Var_epsilon[j, m] <- 1/rgamma(1,
                                    shape = c+ (J/2),
                                    rate = d+ sum( (W[, j] - Mu[,m] - indBeta)^2 )/2 )
      
    }
    
    
    #   # 6. Update W
    
    for (j in 1:J) {
      for (i in 1:N) {
        
        if (r[i, j] == 1) {
          upper <- 100
          lower <- W[which(r[, j] == 2), j]
          
        } else if (r[i,j] == n) {
          upper <- W[which(r[, j] == n - 1), j]
          lower <- -100
          
        } else {
          
          lower <- W[which(r[, j] == r[i,j] + 1), j]
          upper <- W[which(r[, j] == r[i,j] - 1), j]
          
        }
        W[i,j] <- truncnorm::rtruncnorm(1, a = lower, b = upper,
                                        mean = Mu[i, m] + Beta[j]*(i==j),
                                        sd = sqrt( Var_epsilon[j, m] ))
      }
    }
  }
  
  drop = 1:burnin
  post.mean.mus <- apply(Mu[,-drop], 1, mean)
  
  post.mean.betas <- apply(Beta[,-drop], 1, mean)
  post.median.Var_epsilon <- apply(Var_epsilon[,-drop], 1, median)
  
  post.mean.mu_beta <- mean(Mu_beta[,-drop])
  post.median.Var_beta <- median(Var_beta[,-drop])
  
  
  out <- list(
              post_mean_mus = post.mean.mus,
              post_mean_betas = post.mean.betas, 
              post_mean_mubeta = post.mean.mu_beta, 
              post_median_Varbeta = post.median.Var_beta,
              post_median_VarEps = post.median.Var_epsilon,
              post_rank = rank(-post.mean.mus)
  )
  
  return(out)
}  
#-------------------calculations-----------------------
# 

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




#-------------------metrics----------------------------
# 
calculate_top1 <- function(ranked_list, true_list){
  
  order1 <- order(true_list)
  order2 <- order(ranked_list)
  mean(order2[1] %in% order1[1])
}
calculate_top3 <- function(ranked_list, true_list){
  
  order1 <- order(true_list)
  order2 <- order(ranked_list)
  mean(order2[1:3] %in% order1[1:3])
}

calculate_abs_diff <- function(ranked_list, true_list){sum(abs(ranked_list- true_list))}
calculate_ssq_diff <- function(ranked_list, true_list){sum((ranked_list- true_list)^2)}


#-------------------inverse gamma-------------------

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

#------------------Borda scores--------------------
# 02 

myBorda_scores <- function(input, 
                           n,
                           marg=1, # apply to margin of = 1 rows or 2 columns 
                           indSelfEval=TRUE # indicator: use self-evaluation data or not
){
  
  if(!indSelfEval){diag(input)=NA}
  
  mean= apply(input, marg, function(x) mean(x, na.rm = TRUE))
  geomean= apply(input, marg, geo.mean)
  median= apply(input, marg, function(x) median(x, na.rm = TRUE))
  
  no_na_matrix= replace(input, is.na(input), 0)
  l2norm= apply(no_na_matrix, marg, function(x) norm(x,"2"))
  
  rms= apply(input, marg, function(x) sqrt(mean(abs(x)^2, na.rm = TRUE)))
  
  df= data.frame(mean, geomean, median, l2norm, rms)%>%
    melt()%>%   # warning: No id variables; using all as measure variables
    
    group_by(variable)%>%
    mutate(ranking=rank(value))%>%
    mutate(team=factor(paste0("Team",1:n), levels=paste0("Team",1:n)),
           method=ifelse(indSelfEval, "self", "selfrm"))
  
  return(df)
}
#------------------Data simulation--------------------
# 01 

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
    sigma_epsilons <-  (1/rhos^2 - 1)
    
    # Build epsilon matrix column-by-column, one sigma per ranker (column)
    epsilon <- sapply(sigma_epsilons, function(s) rnorm(n, mean = 0, sd = s))
    
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
        sigma_epsilons = sigma_epsilons  # retained for inspection/debugging
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
    df=df
) {
  regenerate <- TRUE  # Flag to control regeneration loop
  result <- NULL
  
  while (regenerate) {
    # Generate data
    mu <- rt(n, df = df) * sqrt(df / (df - 2)) # Global team effects
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
