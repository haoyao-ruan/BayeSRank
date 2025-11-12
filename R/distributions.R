
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


# Function to estimate a and b parameters of an Inverse Gamma distribution

find_shape_rate <- function(mean, sd){ #for Inverse-Gamma
  
  var <- sd^2
  alpha <- (mean^2) / var + 2
  
  # Calculate beta using the formula for mean
  beta <- mean * (alpha - 1)
  
  return(list(shape = alpha, rate = beta))
  
}





















