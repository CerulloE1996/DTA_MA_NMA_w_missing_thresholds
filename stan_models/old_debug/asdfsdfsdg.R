model{
  
  for(i in 1:ns){
    for(j in 1:2){
      
      n[i,j,1] <- N[i,j] 
      p[i,j,1] <- pr[i,j,1] 
      
      for(t in 2:T[i,j]){       # Note (1) (see below)
        n[i,j,t] <- x[i,j,t-1]  		
        p[i,j,t] <-  pr[i,j,t] / pr[i,j,t-1]	 
      }    
    }
    
    for(t in 1:T[i,2]){   # Note (2) 
      q[i,t] <- (( pow(C[i,t], lambda) - 1 ) / lambda )*(1 - equals(lambda, 0))  +  log(C[i,t])*equals(lambda, 0)  	# Box-Cox transformation
    }
    
    for(j in 1:2){
      for(t in 1:T[i,j]){
        
        x[i,j,t] ~ dbin(p[i,j,t], n[i,j,t])  	     # Likelihood                              
        d[i,j,t] <- (mu[i,j] - q[i,t] ) / s[i,j] 
         logit(pr[i,j,t]) <- min(10, max(-10, d[i,j,t]) )   # Note (3)
        
        # xhat[i,j,t] <- p[i,j,t]*n[i,j,t]  		# Fitted values
        # dev[i,j,t] <- 2*(x[i,j,t]*(log(x[i,j,t]) - log(xhat[i,j,t])) + (n[i,j,t] - x[i,j,t])*(log(n[i,j,t] - x[i,j,t])  - log(n[i,j,t] - xhat[i,j,t])))  # Residual deviance contribution
      }
    }
    
    # Distributions of correlated random effects:
    mu[i,1] ~ dnorm(mean[1], prec[1])
    
    mu[i,2] ~ dnorm(cond.mean.mu[i], cond.prec.mu)
    cond.mean.mu[i] <- mean[2] + (rho_mu*sd[2]/sd[1])*(mu[i,1] - mean[1])
    
    for(j in 1:2){
      cond.mean.s[i,j] <- mean[j+2] + (rho_mu_sigma*sd[j+2]/sd[j])*(mu[i,j] - mean[j])
      logs[i,j] ~ dnorm(cond.mean.s[i,j], cond.prec.s[j])I(-5,)
      s[i,j] <- exp(logs[i,j])	
    }
    
    rd[i] <- sum(dev[i,1,1:T[i,1]]) + sum(dev[i,2,1:T[i,2]])	# Residual deviance study i
    
  }
  
  # # Predictive distributions for random effects:
  # mupred[1] ~ dnorm(mean[1], prec[1])
  # cond.mean.mu.pred <- mean[2] + (rho_mu*sd[2]/sd[1])*(mupred[1] - mean[1])
  # mupred[2] ~ dnorm(cond.mean.mu.pred, cond.prec.mu)
  # for(j in 1:2){	
  #   cond.mean.s.pred[j] <- mean[j+2]  + (rho_mu_sigma*sd[j+2]/sd[j])*(mupred[1] - mean[1])  logspred[j] ~ dnorm(cond.mean.s.pred[j], cond.prec.s[j])
  # }
  # 
  # # Priors:
  # lambda ~ dunif(-3,3)
  # for(r in 1:4){
  #   mean[r] ~ dnorm(0, 0.0001)
  #   sd[r] ~ dunif(0,5)
  #   prec[r] <- pow(sd[r], -2)
  # }	
  # rho_mu ~ dunif(-1,1)
  # rho_mu_sigma ~ dunif(-1,1)
  # 
  # cond.var.mu <-  (1- pow(rho_mu,2))*pow(sd[2], 2)
  # cond.prec.mu <- 1/cond.var.mu
  # 
  # for(j in 1:2){
  #   cond.var.s[j]<-  (1- pow(rho_mu_sigma,2))*pow(sd[j + 2], 2)
  #   cond.prec.s[j] <- 1/cond.var.s[j]
  # }
  
  # Summary estimates and predictive intervals. Note (4): 
  for(k in 1:m){
    
    transthres[k] <-  (( pow(thres[k], lambda) - 1 ) / lambda )*(1 - equals(lambda, 0))   +  log(thres[k])*equals(lambda, 0)	
    
    for(j in 1:2){
      logit(summpr[k,j]) <- (mean[j] - transthres[k] ) / exp(mean[j+2])  # Summary estimates 
      logit(predpr[k,j]) <- (mupred[j] - transthres[k] ) / exp(logspred[j])  # Prediction intervals
    }
  }
  
  resdev <- sum(rd[])    # Total residual de