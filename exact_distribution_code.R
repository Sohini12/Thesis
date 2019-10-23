### Codes for finding the expected duration where inputs are data (years and temperature (temp)),q= chosen threshold
### quantiles and max_dur is the maximum length permitted for durations. max_dur can vary from 30 to 365
exact_dist <- function(data,q,max_dur){
  require(forecast)
  n <- dim(data)[1]
  years <- sort(unique(data$years))
  no_of_years <- length(years)
  temp1 <- data$temp
  ord <- rep(0,no_of_years)
  ### Findind the order for every year
  for(i in 1:no_of_years){
    temp <- temp1[which(data$years==years[i])]
    trend <- ma(temp, order=30,centre=T)
    detrend= temp-trend
    ord[i]<- ar(detrend[which(is.na(detrend)==FALSE)], aic = TRUE, order.max = NULL)$order
  }
  
  ### Choosing median of all the orders to be maximum order that we will check
  max_order <- floor(median(ord))
  
  ### Choosing the k with minimum chi-sqaure statistic
  e_array <- array(rep(0,max_dur*1000*max_order),dim=c(1000,max_dur,max_order))
  c <- matrix(rep(0,(max_order*1000)),1000,max_order)  
  med <- rep(0,max_order)
  if(max_order >=4){
    for(o in 4:max_order){
      no_of_order <- o
      freq <- matrix(0,no_of_years,max_dur)
      for(i in 1:no_of_years){
        temp <- temp1[which(data$years==years[i])]
        M_0 <- quantile(temp,q)
        b=c(0,as.numeric(temp>M_0),0)
        d=diff(b)
        t.start=which(d==1)
        t.end=which(d==-1)
        m=min(length(t.start),length(t.end))
        duration=t.end[1:m]-t.start[1:m]
        if(is.na(sum(duration))==TRUE){
          duration=0 }
        
        for(l in 1:max_dur){
          freq[i,l] <- sum(duration==l)
          
        } 
        
      }
      N <- rep(0,no_of_years)
      for(j in 1:no_of_years){
        N[j] <- sum(freq[j,])
        
      }
      
      
      
      
      library('rjags')
      modelstring="
      model{
      
      for(j in 1:no_of_years){
      freq[j,] ~ dmulti(prob[j,],N[j])
      
      }
      
      for(j in 1:no_of_years){
      for(l in 1:(no_of_order-1)){
      prob[j,l] <- pi[j,l]
      
      
      }
      for(l in no_of_order:max_dur){
      prob[j,l] <- A[j]*(B[j]^(l-no_of_order))
      }
      
      }
      
      
      
      
      for(j in 1:no_of_years){
      p[j] ~ dbeta((mu_p*tau_p),((1-mu_p)*tau_p))
      pi[j,1] <- (1-p[j])*alpha[j,1]
      for(i in 2:(no_of_order-2)){
      pi[j,i] <- pi[j,(i-1)]*alpha[j,i]*(1-alpha[j,(i-1)])/alpha[j,(i-1)]
      }
      pi[j,(no_of_order-1)] <- pi[j,(no_of_order-2)]*(1-alpha[j,(no_of_order-2)])/alpha[j,(no_of_order-2)]
      
      }
      for(j in 1:no_of_years){
      for(l in 1:(no_of_order-2)){
      alpha[j,l] ~ dbeta((mu_alpha[l]*tau_alpha[l]),((1-mu_alpha[l])*tau_alpha[l]))T(0.001,0.999)
      
      
      }
      }
      for(l in 1:(no_of_order-2)){
      mu_alpha[l] ~dunif(0.01,0.99)
      tau_alpha[l] ~ dunif(10,40)
      }
      for(j in 1:no_of_years){
      B[j] ~ dbeta((mu_B*tau_B),((1-mu_B)*tau_B))T(0.001,0.999)
      A[j] <- p[j]*(1-B[j])/(1-B[j]^(max_dur+1-no_of_order))
      }
      mu_p ~ dunif(.01,0.99)
      mu_B ~ dunif(0.01,0.99)
      tau_p ~ dunif(10,40)
      tau_B ~ dunif(10,40)
      
      
      
      }
      "
      
      library(rjags)
      jags <- jags.model(textConnection(modelstring),data = list('N'=N,'freq' = freq,'no_of_order'=no_of_order,'no_of_years'=no_of_years,'max_dur'=max_dur))
      
      mu_B_sample <- jags.samples(jags,c('mu_B'),1000)$mu_B
      mu_alpha_sample <- jags.samples(jags,c('mu_alpha'),1000)$mu_alpha
      mu_p_sample <- jags.samples(jags,c('mu_p'),1000)$mu_p
      tau_B_sample <- jags.samples(jags,c('tau_B'),1000)$tau_B
      e <- matrix(rep(0,1000*max_dur),1000,max_dur)
      for(i in 1:1000){
        e[i,1] <- mu_alpha_sample[1,i,]*(1-mu_p_sample[i])
        for(j in 2:(no_of_order-2)){
          e[i,j] <- e[i,(j-1)]*mu_alpha_sample[j,i,]*(1-mu_alpha_sample[(j-1),i,])/mu_alpha_sample[(j-1),i,]
        }
        e[i,(no_of_order-1)] <- e[i,(no_of_order-2)]*(1-mu_alpha_sample[(no_of_order-2),i,])/mu_alpha_sample[(no_of_order-2),i,]
        
      }
      
      
      for(j in no_of_order:max_dur){
        for(i in 1:1000){
          f <- function(b){
            c <- ((1-b)*(b^(j-no_of_order)))*dbeta(b,mu_B_sample[i]*tau_B_sample[i],(1-mu_B_sample[i])*tau_B_sample[i])/(1-b^(max_dur+1-no_of_order))
          }
          e[i,j] <- (mu_p_sample[i])*as.numeric(integrate(f,0,1)[1])
        }
      }
      
      e_array[,,o] <- e
      chi2 <- matrix(rep(0,1000*no_of_years),1000,no_of_years)
      
      for(j in 1:no_of_years){
        for(l in 1:max_dur){
          
          chi2[,j] <- chi2[,j]+ (freq[j,l]-N[j]*e[,l])^2/(N[j]*e[,l] ) 
        }
      }
      
      chi2 <- chi2/(max_dur-no_of_order-1)
      c[,o] <- rowMeans(chi2)
      med[o] <- median(c[,o])
    } 
    
}
  
  
  if(max_order >=2){
    o <- 2
    no_of_order <- o
    
    freq <- matrix(0,no_of_years,max_dur)
    
    for(i in 1:no_of_years){
      temp <- temp1[which(data$years==years[i])]
      M_0 <- quantile(temp,q)
      b=c(0,as.numeric(temp>M_0),0)
      d=diff(b)
      t.start=which(d==1)
      t.end=which(d==-1)
      m=min(length(t.start),length(t.end))
      duration=t.end[1:m]-t.start[1:m]
      if(is.na(sum(duration))==TRUE){
        duration=0 }
      
      for(l in 1:max_dur){
        freq[i,l] <- sum(duration==l)
        
      } 
      
    }
    N <- rep(0,no_of_years)
    for(j in 1:no_of_years){
      N[j] <- sum(freq[j,])
      
    }
    
    
    
    
    library('rjags')
    modelstring="
    model{
    
    for(j in 1:no_of_years){
    freq[j,] ~ dmulti(prob[j,],N[j])
    
    }
    
    for(j in 1:no_of_years){
    for(l in 1:(no_of_order-1)){
    prob[j,l] <- pi[j,l]
    
    
    }
    for(l in no_of_order:max_dur){
    prob[j,l] <- A[j]*(B[j]^(l-no_of_order))
    }
    
    }
    
    
    
    
    for(j in 1:no_of_years){
    p[j] ~ dbeta((mu_p*tau_p),((1-mu_p)*tau_p))T(0.001,0.999)
    pi[j,1] <- (1-p[j])
    
    
    
    }
    
    
    
    
    for(j in 1:no_of_years){
    B[j] ~ dbeta((mu_B*tau_B),((1-mu_B)*tau_B))T(0.001,0.999)
    A[j] <- p[j]*(1-B[j])/(1-B[j]^(max_dur+1-no_of_order))
    }
    mu_p ~ dunif(0.01,0.99)
    mu_B ~ dunif(0.01,0.99)
    tau_p ~ dunif(10,40)
    tau_B ~ dunif(10,40)
    
    
    
    }
    "
    
    library(rjags)
    jags <- jags.model(textConnection(modelstring),data = list('N'=N,'freq' = freq,'no_of_order'=no_of_order,'no_of_years'=no_of_years,'max_dur'=max_dur))
    
    mu_B_sample <- jags.samples(jags,c('mu_B'),1000)$mu_B
    
    mu_p_sample <- jags.samples(jags,c('mu_p'),1000)$mu_p
    tau_B_sample <- jags.samples(jags,c('tau_B'),1000)$tau_B
    e <- matrix(rep(0,1000*max_dur),1000,max_dur)
    for(i in 1:1000){
      e[i,1] <- (1-mu_p_sample[i])
    }
    
    for(j in no_of_order:max_dur){
      for(i in 1:1000){
        f <- function(b){
          c <- ((1-b)*(b^(j-no_of_order)))*dbeta(b,mu_B_sample[i]*tau_B_sample[i],(1-mu_B_sample[i])*tau_B_sample[i])/(1-b^(max_dur+1-no_of_order))
        }
        e[i,j] <- (mu_p_sample[i])*as.numeric(integrate(f,0,1)[1])
      }
    }
    
    e_array[,,o] <- e
    chi2 <- matrix(rep(0,1000*no_of_years),1000,no_of_years)
    
    for(j in 1:no_of_years){
      for(l in 1:max_dur){
        for(i in 1:1000){
          chi2[i,j] <- chi2[i,j]+ (freq[j,l]-N[j]*e[i,l])^2/(N[j]*e[i,l] ) 
        }
      }
    }
    chi2 <- chi2/(max_dur-no_of_order-1)
    c[,o] <- rowMeans(chi2)
    med[o] <- median(c[,o])
  }
  
  if(max_order >=3){
    o <- 3
    no_of_order <- o
    
    freq <- matrix(0,no_of_years,max_dur)
    
    for(i in 1:no_of_years){
      temp <- temp1[which(data$years==years[i])]
      M_0 <- quantile(temp,q)
      b=c(0,as.numeric(temp>M_0),0)
      d=diff(b)
      t.start=which(d==1)
      t.end=which(d==-1)
      m=min(length(t.start),length(t.end))
      duration=t.end[1:m]-t.start[1:m]
      if(is.na(sum(duration))==TRUE){
        duration=0 }
      
      for(l in 1:max_dur){
        freq[i,l] <- sum(duration==l)
        
      } 
      
    }
    N <- rep(0,no_of_years)
    for(j in 1:no_of_years){
      N[j] <- sum(freq[j,])
      
    }
    
    
    
    
    library('rjags')
    modelstring="
    model{
    
    for(j in 1:no_of_years){
    freq[j,] ~ dmulti(prob[j,],N[j])
    
    }
    
    for(j in 1:no_of_years){
    for(l in 1:(no_of_order-1)){
    prob[j,l] <- pi[j,l]
    
    
    }
    for(l in no_of_order:max_dur){
    prob[j,l] <- A[j]*(B[j]^(l-no_of_order))
    }
    
    }
    
    
    
    
    for(j in 1:no_of_years){
    p[j] ~ dbeta((mu_p*tau_p),((1-mu_p)*tau_p))
    pi[j,1] <- (1-p[j])*alpha[j,1]
    
    pi[j,(no_of_order-1)] <- (1-p[j])*(1-alpha[j,1])
    
    }
    for(j in 1:no_of_years){
    for(l in 1:(no_of_order-2)){
    alpha[j,l] ~ dbeta((mu_alpha[l]*tau_alpha[l]),((1-mu_alpha[l])*tau_alpha[l]))T(0.001,0.999)
    
    
    }
    }
    for(l in 1:(no_of_order-2)){
    mu_alpha[l] ~dunif(0.01,0.99)
    tau_alpha[l] ~ dunif(10,40)
    }
    for(j in 1:no_of_years){
    B[j] ~ dbeta((mu_B*tau_B),((1-mu_B)*tau_B))T(0.001,0.999)
    A[j] <- p[j]*(1-B[j])/(1-B[j]^(max_dur+1-no_of_order))
    }
    mu_p ~ dunif(0.01,0.99)
    mu_B ~ dunif(0.01,0.99)
    tau_p ~ dunif(10,40)
    tau_B ~ dunif(10,40)
    
    
    
    }
    "
    
    library(rjags)
    jags <- jags.model(textConnection(modelstring),data = list('N'=N,'freq' = freq,'no_of_order'=no_of_order,'no_of_years'=no_of_years,'max_dur'=max_dur))
    
    mu_B_sample <- jags.samples(jags,c('mu_B'),1000)$mu_B
    mu_alpha_sample <- jags.samples(jags,c('mu_alpha'),1000)$mu_alpha
    mu_p_sample <- jags.samples(jags,c('mu_p'),1000)$mu_p
    tau_B_sample <- jags.samples(jags,c('tau_B'),1000)$tau_B
    e <- matrix(rep(0,1000*max_dur),1000,max_dur)
    for(i in 1:1000){
      e[i,1] <- mu_alpha_sample[1,i,]*(1-mu_p_sample[i])
      
      e[i,(no_of_order-1)] <- (1-mu_alpha_sample[1,i,])*(1-mu_p_sample[i])
      
    }
    for(j in no_of_order:max_dur){
      for(i in 1:1000){
        f <- function(b){
          c <- ((1-b)*(b^(j-no_of_order)))*dbeta(b,mu_B_sample[i]*tau_B_sample[i],(1-mu_B_sample[i])*tau_B_sample[i])/(1-b^(max_dur+1-no_of_order))
        }
        e[i,j] <- (mu_p_sample[i])*as.numeric(integrate(f,0,1)[1])
      }
    }
    
    e_array[,,o] <- e
    chi2 <- matrix(rep(0,1000*no_of_years),1000,no_of_years)
    
    for(j in 1:no_of_years){
      for(l in 1:max_dur){
        
        chi2[,j] <- chi2[,j]+ (freq[j,l]-N[j]*e[,l])^2/(N[j]*e[,l] ) 
      }
    }
    
    chi2 <- chi2/(max_dur-no_of_order-1)
    c[,o] <- rowMeans(chi2)
    med[o] <- median(c[,o])
  }
  
  min_k <- which.min(med[2:max_order])+1
  e <- e_array[,,min_k]
  
  ### Computing the posterior distribution of expected duration
  g <- matrix(rep(0,1000*max_dur),1000,max_dur)
  for(l in 1:max_dur){
    g[,l] <- l*e[,l]
  }
  exp_d_dist<- rowSums(g)
  ### Computing posterior median and credible interval
  exp_d<- median(exp_d_dist)
  
  
  up_est_exp_d <- quantile(exp_d_dist,0.975)
  low_est_exp_d <- quantile(exp_d_dist,0.025)
  ### Output
  fit <- NULL
  fit$no_of_years <- no_of_years  
  ### expected duration
  fit$duration <- exp_d
  
  ### Posterior samples of expected durations
  fit$duration_posterior <- exp_d_dist
  
  ### Credible intervals
  fit$upper <- up_est_exp_d
  fit$lower <- low_est_exp_d
  ### Order that is finally fitted
  fit$order <- min_k
  fit$expected <- e
  return(fit)
  }