# ──────────────────────────────────────────────────────────────────────────────
# Two–step BayesRank 
# ──────────────────────────────────────────────────────────────────────────────
bayesrank_2step <- function(rank_matrix,
                            M_itrns      = 2000,
                            m_burn       = 1000,
                            init_W       = NULL,
                            # ↓ optional pilot‑stage hyper‑priors
                            step1_minbeta      = -100,
                            step1_maxbeta      =  100,
                            step1_Var_beta0    = 1,
                            step1_Var_epsilon0 = 1,
                            step1_a = 0.1, step1_b = 0.1,
                            step1_c = 0.1, step1_d = 0.1,
                            # ↓ random‑seed helper so runs are reproducible
                            seed = 1111) {
  
  if (!is.null(seed)) set.seed(seed)
  
  ## --------------------------------------------------------------------------
  ## 0.  Build or accept an initial W  (un‑ordered latent scores)
  ## --------------------------------------------------------------------------
  n_teams <- nrow(rank_matrix)
  
  if (is.null(init_W)) {
    init_W <- matrix(rnorm(n_teams^2, mean = 0, sd = 1), nrow = n_teams)
  }
  
  ## Re–order each column of W so that larger latent utility ⇒ higher rank
  sorted_W <- apply(init_W, 2, sort, decreasing = TRUE)
  initW <- mapply(function(sorted_col, rank_col) sorted_col[rank_col],
                  as.data.frame(sorted_W),
                  as.data.frame(rank_matrix))
  initW <- matrix(unlist(initW), nrow = n_teams)
  
  ## --------------------------------------------------------------------------
  ## 1.  Preliminary chain with weak / diffuse priors
  ## --------------------------------------------------------------------------
  prelim <- BayeSRank(n = n_teams,
                      r = rank_matrix,
                      M = M_itrns,
                      burnin = m_burn,
                      W = initW,
                      mu0         = rnorm(n_teams, 0, 1),
                      mubeta0     = rnorm(1, 0, 1),
                      Var_beta0   = step1_Var_beta0,
                      Var_epsilon0= step1_Var_epsilon0,
                      minbeta     = step1_minbeta,
                      maxbeta     = step1_maxbeta,
                      a = step1_a, 
                      b = step1_b,
                      c = step1_c, 
                      d = step1_d)
  
  ## Empirical tuning from pilot output
  Beta.bar <- mean(prelim$Beta)
  Beta.sd  <-  sd(prelim$Beta)
  minBeta_est <- Beta.bar - 3 * Beta.sd
  maxBeta_est <- Beta.bar + 3 * Beta.sd
  varBeta_est    <- median(prelim$Var_beta)
  varEpsilon_est <- median(prelim$Var_epsilon)
  
  IG_Beta_est    <- find_shape_rate(varBeta_est,    sd(prelim$Var_beta))
  IG_Epsilon_est <- find_shape_rate(varEpsilon_est, sd(prelim$Var_epsilon))
  
  
  ## --------------------------------------------------------------------------
  ## 2.  Main chain with tuned hyper‑priors
  ## --------------------------------------------------------------------------
  res <- BayeSRank(n = n_teams,
                   r = rank_matrix,
                   M = M_itrns,
                   burnin = m_burn,
                   W = initW,
                   mu0         = rnorm(n_teams, 0, 1),
                   mubeta0     = rnorm(1, 0, 1),
                   Var_beta0   = varBeta_est,
                   Var_epsilon0= varEpsilon_est,
                   minbeta     = minBeta_est,
                   maxbeta     = maxBeta_est,
                   a = IG_Beta_est$shape, 
                   b = IG_Beta_est$rate,
                   c = IG_Epsilon_est$shape, 
                   d = IG_Epsilon_est$rate)
  
  
    return(res)
}


## --------------------------------------------------------------------------
## BayeSRank function
## --------------------------------------------------------------------------
BayeSRank <- function(n, 
                      r,
                      M,
                      burnin, 
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
    
    
    #   # 2.Update Beta
    
    for (j in 1:J) { #for each ranker in COLUMNS
      
      g <- ( (W[j, j] - Mu[j, m]) / Var_epsilon[j, m-1] ) + (Mu_beta[m-1]/Var_beta[m-1]) #numerator g
      h <- (1 / Var_epsilon[j, m-1]) + (1/Var_beta[m-1])              #denominator h
      
      Beta[j, m] <- truncnorm::rtruncnorm(1, a = minbeta, b = maxbeta,
                                          mean = (g/h),
                                          sd = sqrt(1/h))
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
          upper <- Inf
          lower <- W[which(r[, j] == 2), j]
          
        } else if (r[i,j] == n) {
          upper <- W[which(r[, j] == n - 1), j]
          lower <- -Inf
          
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
  cor <- cor(rank(-post.mean.mus), true_rank, use = "everything", method="spearman")
  
  post.mean.betas <- apply(Beta[,-drop], 1, mean)
  
  post.mean.mu_beta <- mean(Mu_beta[,-drop])
  post.median.Var_beta <- median(Var_beta[,-drop])
  
  top1rate <- calculate_top1(ranked_list=rank(-post.mean.mus), true_list=true_rank)
  top3rate <- calculate_top3(ranked_list=rank(-post.mean.mus), true_list=true_rank)
  
  out <- list(# Mu = Mu, 
    Beta = Beta[,-drop],
    #Mu_beta = Mu_beta[,-drop],
    Var_beta = Var_beta[,-drop],
    Var_epsilon = Var_epsilon[,-drop],
    # W=W,
    post_mean_mus = post.mean.mus,
    # post_mean_betas = post.mean.betas, 
    post_mean_mubeta = post.mean.mu_beta, 
    post_median_Varbeta = post.median.Var_beta,
    corr_spearman = cor,
    top1rate=top1rate,
    top3rate=top3rate
    # post_rank = rank(-post.mean.mu),
    # true_rank = true_rank,
    # true_beta = mydat$true_beta
  )
  
  return(out)
}  


