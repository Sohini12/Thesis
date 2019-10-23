### Code for finding the expected duration using the approximate distribution where inputs 
### are data (years and temperature (temp)) and quant is the chosen threshold quantiles.
Duration_func <- function(data,quant){
  
  
  temp1 <- data$temp
  years <- sort(unique(data$year))
  no_of_years <- length(years)
  M_0 <- quantile(temp1,quant)
  b1=c(0,as.numeric(temp1>M_0),0)
  
  mle_estimate <- rep(0,no_of_years)
  no_of_exceedances <- rep(0,no_of_years)
  
  for(i in 1:no_of_years){
    temp <- temp1[grep(years[i],data$year)]
    b=c(0,as.numeric(temp>M_0),0)
    d=diff(b)
    t.start=which(d==1)
    t.end=which(d==-1)
    m=min(length(t.start),length(t.end))
    duration=t.end[1:m]-t.start[1:m]
    no_of_exceedances[i] <- length(duration)
    mle_estimate[i] <- 1-(1/mean(duration))
    
  }
  
  maximum <- max(no_of_exceedances)
  duration <- matrix(rep(0,(maximum*no_of_years)),maximum,no_of_years)
  
  
  for(i in 1:no_of_years){
    temp <- temp1[grep(years[i],data$year)]
    b=c(0,as.numeric(temp>M_0),0)
    d=diff(b)
    t.start=which(d==1)
    t.end=which(d==-1)
    m=min(length(t.start),length(t.end))
    duration[1:m,i]=t.end[1:m]-t.start[1:m]
  }
  
  LL <- function(alpha,tau) {
    R2 <- 0
    for(j in 1:dim(duration)[2]){
      R2 <- R2+lbeta((sum(duration[,j][which(duration[,j]!=0)])-length(duration[,j][which(duration[,j]!=0)])+alpha*tau),(length(duration[,j][which(duration[,j]!=0)])+(1-alpha)*tau))
    }
    
    -R2+no_of_years*(lbeta(alpha*tau,(1-alpha)*tau))
  }
  
  library(stats4)
  
  options(warn=-1)
  a <-mle(LL,start = list(alpha = 0.5, tau = 1))
  alpha_est <- coef(summary(a))[1,1]
  tau_est <- coef(summary(a))[2,1]
  variance_covariance_est <- vcov(a)
  
  
  p_value_ks <- ks.test(mle_estimate[which(is.na(mle_estimate)==FALSE)], "pbeta",(alpha_est*tau_est), (1-alpha_est)*tau_est)$p
  require(goftest)
  p_value_ad <- ad.test(mle_estimate[which(is.na(mle_estimate)==FALSE)], "pbeta",(alpha_est*tau_est), (1-alpha_est)*tau_est)$p
  
  
  fit <- NULL
  fit$no_of_years <- no_of_years
  fit$alpha <- alpha_est
  fit$tau <- tau_est
  fit$var_covar <- variance_covariance_est
  fit$p_value_ks <- p_value_ks
  fit$p_value_ad <- p_value_ad
  fit$duration <- mle_estimate[which(is.na(mle_estimate)==FALSE)]
  
  return(fit)
}

