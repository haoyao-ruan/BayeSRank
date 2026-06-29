# ============================================================================
# competing_methods.R — baseline / competing aggregators
# Borda-score family, BIRRA and BiGER.
# Part of the BayeSPeer function library; logic unchanged from the original implementation.
# ============================================================================

## --- Borda-score aggregators (Mean / Geomean / Median / L2norm) ------------
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

## --- BIRRA -----------------------------------------------------------------
BIRRA=function(data, prior=0.05, num.bin=50, num.iter=10, return.all=F, plot.legend=F, grp=NULL, cor.stop=1, ...){
  nr=nrow(data)
  nrp=floor(nrow(data)*prior)
  #data=apply(-data,2,rank)/nr
  data=apply(data,2,rank)/nr
  
  nc=ncol(data)
  TPR=FPR=Bayes.factors=matrix(ncol=nc, nrow=num.bin)
  binned.data=ceiling(data*num.bin)
  
  
  bayes.data=matrix(nrow=nrow(data), ncol=ncol(data))
  
  guess=apply(data,1,mean)
  cprev=0
  #par(mfrow=c(floor(sqrt(num.iter)), ceiling(sqrt(num.iter))), mai=rep(0.7,4))
  for ( iter in 1:num.iter){
   # if((cor.stop-cprev)>1e-15){
     if((cor.stop-cprev)>1e-15| is.na(cprev)){ #updated by Haoyao for selfrm data
       
      guesslast=guess
      oo=order(guess)
      guess[oo[1:nrp]]=1
      guess[oo[(nrp+1):nr]]=0
      
      
      for (i in 1:nc){  
        for (bin in 1:num.bin){
          frac=bin/num.bin
          TPR=sum(guess[binned.data[,i]<=bin])
          FPR=sum((!guess)[binned.data[,i]<=bin])
          
          Bayes.factors[bin,i]=log((TPR+1)/(FPR+1)/(prior/(1-prior)))
          
        }
      }
      
      Bayes.factors=apply(Bayes.factors,2,smooth)
      Bayes.factors=apply(Bayes.factors,2,function(x){rev(cummax(rev(x)))})
      # Plot TPR vs bin for each data set
      if(is.null(grp)){
        matplot(1:num.bin, Bayes.factors, type="l", lwd=2, ...)
      }
      else{
        matplot(1:num.bin, Bayes.factors, type="l", lwd=2, lty=grp, col=grp)
      }
      
      title(paste("Iteration", iter))
      if (iter==1&plot.legend){
        legend("topright", col=1:5, lty=1:4, legend=colnames(data), lwd=2, ncol=2)
      }
      for (bin in 1:num.bin){
        oo=order(Bayes.factors[bin,], decreasing=T)
        Bayes.factors[bin, oo[1]]=Bayes.factors[bin, oo[2]]
        
        
      }
      
      for (i in 1:nc){
        
        bayes.data[,i]=Bayes.factors[binned.data[,i],i]
        
      }
      
      
      bb=exp(apply(bayes.data,1, sum))
      f=prior/(1-prior)
      prob=bb*f/(1+bb*f)
      exp=sort(prob, decreasing=F)[nrp]
      
      guess=rank(-apply(bayes.data,1, sum))
      cprev=cor(guess, guesslast)
      message("correlation with pervious iteration=",cprev)
    }
    else{
      message("Converged");
      break
    }
  }
  if(return.all){
    return(list(result=guess, data=bayes.data, BF=Bayes.factors))
  }
  else{
    guess
  }
}

## --- BiGER -----------------------------------------------------------------

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