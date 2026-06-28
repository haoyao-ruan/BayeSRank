# ============================================================================
# sampler.R — BayeSPeer Gibbs samplers (Inverse-Gamma priors)
# Pilot (.pre) and main Gibbs samplers.
# Part of the BayeSPeer function library. Functions are moved verbatim from the
# original 00_functions.R / method scripts; logic is unchanged.
# ============================================================================

#------------------Gibbs Sampler Application-------------------
BayeSPeer.pre <- function(n, r, M, burnin, 
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
BayeSPeer <- function(n, r, M, burnin, 
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
