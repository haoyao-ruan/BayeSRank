
## BiGER Algorithm
BiGER <- function(r,
                  W,
                  Beta0=rnorm(ncol(r), 0, 1),
                  tau2=100,
                  sigma_s0=rep(1,ncol(r)),
                  sigma_sb0 = 1,
                  a=.1,
                  b=.1,
                  c=.1,
                  d=.1,
                  mu0=numeric(nrow(r)),
                  M = 5000, burnin = 2500,
                  verbose=1000){
  # Arrays
  I <- dim(r)[1]
  J <- dim(r)[2]
  Mu <- matrix(0,I,M)
  Sigma_s2 <- matrix(0,J,M)
  Sigma_sb <- rep(0,M)
  Mu_beta <- rep(0, M)
  Beta <- matrix(0, J, M)
  
  # Inits
  Beta[,1] <- Beta0
  Mu[,1] <- mu0
  Sigma_s2[,1] <- sigma_s0
  Sigma_sb[1] <- sigma_sb0
  
  for (m in 2:M)
  {
    if (m %% verbose == 0) {
      print(m)
    } 
    for (i in 1:I)
    {
      b_temp <- rep(0, J)
      b_temp[i] <- Beta[i, m-1]
      numerator <- sum((W[i,]-b_temp)/Sigma_s2[,m-1])
      denominator <- sum(1/Sigma_s2[,m-1])+1
      Mu[i,m] <- rnorm(1,numerator/denominator,sqrt(1/d))
    }
    
    # Update Mu Beta
    numerator <- sum(Beta[,m-1])/Sigma_sb[m-1]
    denominator <- 1/Sigma_sb[m-1] + 1/tau2
    Mu_beta[m] <- rnorm(1,numerator/denominator,sqrt(1/d))
    
    # Update Beta
    for (j in 1:J) {
      numerator <- (W[i,i]-Mu[i, m])/Sigma_s2[i,m-1]
      denominator <- 1/Sigma_s2[i,m-1]+ 1/Sigma_sb[m-1]
      Beta[j,m] <- rnorm(1,numerator/denominator,sqrt(1/d))
    }
    
    ## updata sigma_s2
    for (j in 1:J)
    {
      Sigma_s2[j,m]=1/rgamma(1,shape=J/2+a,
                             rate=sum((W[,j]-Mu[,m])^2)/2+b)
    }
    
    ## updata sigma_sb
    Sigma_sb[m]=1/rgamma(1,shape=J/2+c,
                         rate=sum((Beta[,m]-Mu_beta[m])^2)/2+d)
    
    #update W
    for (j in 1:J) {
      sigma_s2_sm <- sqrt(Sigma_s2[j,m])
      for (i in 1:I) {
        if (r[i, j]==1) {
          lower = W[which(r[,j]==2),j]
          upper = Inf
        } else if (r[i, j]==I) {
          lower=-Inf
          upper = W[which(r[,j]==(I-1)),j]
        } else {
          lower=W[which(r[,j]==r[i,j]+1),j]
          upper = W[which(r[,j]==r[i,j]-1),j]
        }
        
        W[i,j]=truncnorm::rtruncnorm(1,a=lower,b=upper,
                                     mean=Mu[which(r[,j]==i),m]+Beta[i, m]*(i==j),
                                     sd=sigma_s2_sm) 
        
      }
    }
  }
  post.mean.mu=base::apply(Mu[,ceiling(burnin):M],1,mean)
  
  out <- list(post.mean.mu=post.mean.mu, Sigma_s2 =Sigma_s2 , Sigma_s=Sigma_sb)
  return(out)
}