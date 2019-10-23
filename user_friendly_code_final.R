#### General functions

j.fun.tau=function(tau,resids,sigma.2){
  # Function to compute the unpenalized function for finding the saddle point
  # Args:
  #     tau: tau in the function
  #  resids: residuals obtained using corresponding beta
  # sigma.2: square-root of the target residual variance (can be obtained from 
  #                                                 initial robust estimation)

  mean(exp(-tau*(resids^2-sigma.2)))
} 



ETReg_new<-function(x,y,lambda,beta_previous,tau.hat){
  # Function to choose the estimate of beta
  #  Args:
  #             x: the desigm matrix with all the predictors (without the intercept)
  #             y: the vector of response variable
  #        lambda: a single choice of lambda
  #       tau.hat: an estimate of tau
  # beta_previous: an initial choice of beta
  require(MASS)
  require(glmnet)
  require(matrixcalc)
  n=length(y)   #number of observations
  p = dim(x)[2] #number of predictors
  
  lts.fits=sparseLTS(x,y)
  lasso.fit <- glmnet(x, y,lambda= lambda,intercept = TRUE)
  beta.01 <- matrix(coef(lasso.fit)) 
  # choice of another initial beta by LASSO regression
  sigma.2=lts.fits$scale^2
  # target residual variance 
  hat.beta=cbind(beta_previous,beta.01)
  hat.tau=NULL
  
  beta=hat.beta +1
  
  l1_norm <- rep(0,2)
 
 
  iter=2
  beta_conv <-  matrix(0,10001,p+1)
  
  for (start.val in 1:2){
    iter=2
    beta_conv <-  matrix(0,1001,p+1)
    
    while (sum(abs(beta-hat.beta[,start.val]))/sum(abs(beta))>0.000001 && iter<1001)
    {
      
      beta = hat.beta[,start.val]
      for(pred in 1:p){
        x1=cbind(rep(1,n),x)
        resids=y-x1%*%beta
        single_deriv <- matrix(0,(n+1),1)
        single_deriv[1] <- 0
        for(k in 2:(n+1)){
          single_deriv[k,] <-  single_deriv[k-1]+ exp(-tau.hat*(resids[k-1]^2-sigma.2))*2*tau.hat*resids[k-1]*x[k-1,pred]
        }
        single_derivative <- single_deriv[n+1]
        double_deriv <- matrix(0,(n+1),1)
        double_deriv[1] <- 0
        for(k2 in 2:(n+1)){
          double_deriv[k2] <- double_deriv[k2-1]+exp(-tau.hat*(resids[k2-1]^2-sigma.2))*(4*(tau.hat^2)*resids[k2-1]^2- 2*tau.hat)*x[k2-1,pred]^2
        }
        double_derivative <- double_deriv[n+1]
        if(double_derivative < 0){
          v <- sqrt((-1/2)*double_derivative)
          u <- (single_derivative - double_derivative*beta[pred+1])/v
          b_ols <- lm(u~v-1)$coef
          l <- lambda/(2*v^2)
          if (b_ols > l){
            beta[pred+1] <- b_ols - l
          }
          if(b_ols < -l){
            beta[pred+1] <- b_ols +l
          }
          if(b_ols <=l && b_ols>= -l){
            beta[pred+1] <- 0
          } 
        }
        else{
          
          stepsize <- 1
          beta[pred+1] <- beta[pred+1]+ stepsize*single_derivative
          
        }
        hat.beta[,start.val] <- beta
        beta_conv[iter,] <- beta
        iter = iter +1
      }
    }
    l1_norm[start.val] <- sum(abs(hat.beta[,start.val]))
    
  }
    
  final_beta <- hat.beta[,which.min(l1_norm)]
  fit=NULL
  fit$coef <- final_beta
  return(fit)
}





#### Simulation Codes for proposed method
### Simulating data
## 1. with no outliers

set.seed(seed) 
# Assigning a seed before every simulation
n <- 100
p <- 40
z <- mvrnorm(n, rep(0,p),diag(p))
z1 <- rnorm(n)
x <- matrix(0,n,p)
for(i in 1:n){
  for (j in 1:p){
    x[i,j] <- z[i,j] + z1[i]
  }}
# simulation of predictor variables
# n: number of observations to be simulated (we used n=100)
# p: number of predictors to be simulated (we used p=40)


beta <- c(rep(0,10), rep(2,10),rep(0,10),rep(2,10))
# True beta for our case (can be any vector of length p 
#           with some entries to be exactly equal to 0)

y <- x%*%beta +15*rnorm(n) 
# Simulation of response variables

## 2. with 10 percent vertical outliers

set.seed(seed)
z <- mvrnorm(n, rep(0,p),diag(p))
z1 <- rnorm(n)
x <- matrix(0,n,p)
for(i in 1:n){
  for (j in 1:p){
    x[i,j] <- z[i,j] + z1[i]
  }}
beta <- c(rep(0,10), rep(2,10),rep(0,10),rep(2,10))
y1 <- x[1:(n*0.9),]%*%beta +15*rnorm(n*0.9) 
y2 <- x[(n*0.9+1):n,]%*%beta+15*rnorm((n*0.1),mean= 50)
# 10 % of the outliers in the response variable
y <- c(y1,y2)


## 3. with 10 5 leverage points

set.seed(seed)
z_1 <- mvrnorm((n*0.9), rep(0,p),diag(p))
z_2 <- mvrnorm((n*0.1), rep(5,p),diag(p))
# 10 5 of the outliers in the predictor variables
z <- rbind(z_1,z_2)
z1 <- rnorm(100)
x <- matrix(0,100,40)
for(i in 1:100){
  for (j in 1:40){
    x[i,j] <- z[i,j] + z1[i]
  }}
beta <- c(rep(0,10), rep(2,10),rep(0,10),rep(2,10))
y1 <- x[1:(n*0.9),]%*%beta +15*rnorm(n*0.9) 
y2 <- x[(n*0.9+1):n,]%*%beta+15*rnorm((n*0.1),mean= 50)
# 10 % of the outliers in the response variable
y <- c(y1,y2)


## 4. with 10 % outliers with alternative beta

z <- mvrnorm(n, rep(0,p),diag(p))
z1 <- rnorm(n)
x <- matrix(0,n,p)
for(i in 1:n){
  for (j in 1:p){
    x[i,j] <- z[i,j] + z1[i]
  }}
beta <- c(rep(0,10), rep(2,10),rep(0,10),rep(2,10))
beta1 <- c(rep(7,10),rep(2,10),rep(0,10),rep(2,10))
# Alternative choice of beta
y1 <- x[1:(n*0.9),]%*%beta +15*rnorm((n*0.9)) 
y2 <- x[(n*0.9+1):n,]%*%beta1+3*rnorm((n*0.1))
# 10 % of the outliers from the alternative beta
y <- c(y1,y2)


### Finding the solution path for the proposed method using cross validation

ETReg_sol_path =  function(x,y,beta_intial,scale_initial){
  # Note: this function does not give stepwise predictor variable selection
  #       (more than one lambda can have same number of predictors in the estimate).
  #       To get a stepwise predictor selection, use function ETReg_PM_sol_path 
  #       (described below)
  # Args:
  #   beta_initial: an initial robust estimate of beta
  #  scale_initial: a target residual variance obtained from a robust estimation
  #              x: the design matrix with predictors
  #              y: the vector of response variables
  #       no_of_cv: number of cross-validations to perform
  require(MASS)
  require(glmnet)
  require(matrixcalc)
  n=length(y)
  p = dim(x)[2]
  x = matrix(x, nrow=n)
  beta.01 <- beta_initial
  sigma.2 <-scale_initial^2
  x1=cbind(rep(1,n),x)
  resids=y-x1%*%beta.01
  tau.hat=optimize(j.fun.tau,c(0,1000/sigma.2),resids=resids,sigma.2=sigma.2)$minimum
  # finding the minimizer of tau
  beta_previous <- beta.01
  lambda_grid <- function(x,y,lambda_min,lambda_max,tau.hat,sigma.2){
  # function to find the sequence of lambda to obtain a solution path
  # Args:
  #             x: design matrix of predictors
  #             y: vector of response variable
  #    lambda_min: minimum value of lambda in the sequence
  #    lambda_max: maximum value of lambda in the sequence
  #       sigma.2: sqaure-root of the target residual obtained by intial robust estimation
  #       tau.hat: an estimate of tau
    
    fit_min <- ETReg_new(x,y,lambda_min,beta_initial,tau.hat)
    
    p_min <- sum(fit_min$coef[2:(p+1)]!=0)
    fit_max <- ETReg_new(x,y,lambda_max,fit_min$coef,tau.hat)
    
    p_max <- sum(fit_max$coef[2:(p+1)]!=0)
    
    lambda <- (lambda_min+lambda_max)/2
    fit <- ETReg_new(x,y,lambda,fit_min$coef,tau.hat,sigma.2)
    
    p_fit <- sum(fit$coef[2:(p+1)]!=0)
    
    beta_result <- c(lambda,fit$coef)
    
    
    if( lambda_max- lambda_min< 0.000001){
      
      
      return(beta_result)
      
    }
    else{ 
      if (p_min- p_fit == 0 && p_fit- p_max == 0){
        
        
        return(beta_result)
      }else{     
        
        l <- cbind(lambda_grid(x,y,lambda_min, lambda,beta_previous),lambda_grid(x,y,lambda, lambda_max,fit_min$coef))
        return(l)
      }
    }
    
    
  }
  
  lambda_grid_fit <- lambda_grid(x,y,lambda_min=0,lambda_max=1,tau.hat=tau.hat,sigma.2=sigma.2) 
  beta_final<- t(lambda_grid_fit[1:(p+2),])
  # lambda grid with their estimate of beta
  
  beta_final <- beta_final[order(beta_final[,1]),]
  lambda_series <- beta_final[,1]
  #lambda grid
  beta_final <- beta_final[,-1]
  beta_final <- as.matrix(beta_final)
  #solution path of beta
  fit <- NULL
  fit$sol_path <- beta_final
  fit$lambda_grid <- lambda_series
  return(fit)
}


ETReg_solution = function(x,y,k,h,beta_initial,scale_initial,seed){
  # function to find the optimal lambda and corresponding solution 
  #                                         using cross-validation
  #        Args:
  #             x: desigm matrix with predictor variables
  #             y: vector of response variable
  #             k: number of folds of cross-validation
  #             h: trimming proportion
  #  beta_initial: an initial choice of beta
  # scale_initial: square-root of the target residual variance obtained
  #                                    by an initial robust estimation
  #         seed: seed assignment for cross-validation
  fit_sol_path <- ETReg_sol_path(x,y,beta_initial,scale_initial)
  beta_final <- fit_sol_path$sol_path
  length <- dim(beta_final)[1]
  res <- matrix(0,no_of_cv,length)
  lambda_series <- fit_sol_path$lambda_grid
  n <- length(y)
  p <- dim(x)[2]
  for(j in 1:no_of_cv){
    split <- sample(c(1:n),n)
    x_split <- x[split,]
    y_split <- y[split]
    test <- array(0,c(k,(n/k),p))
    test_y <- array(0,c(k,1,(n/k)))
    train <- array(0,c(k,(n-(n/k)),p))
    train_y <- array(0,c(k,1,(n-(n/k))))
    beta_train <- array(0,c(k,(p+1),1))
    residual <- array(0,c(k,1,(n/k)))
    tau.hat <- rep(0,k)
   
    for(l in 1:length){
      
      lambda <- lambda_series[l]
      
      for(i in 1:k){
        
       
        beta.01 <- beta_initial
        sigma_2 <-scale_initial
        x1=cbind(rep(1,(n-(n/k))),train[i,,])
        resids=y-x1%*%beta.01
        tau.hat[i]=optimize(j.fun.tau,c(0,1000/sigma_2),resids=resids,sigma.2=sigma_2)$minimum
        test[i,,] <- x_split[(1+(i-1)*(n/k)):(i*(n/k)),]
        
        test_y[i,,] <- y_split[(1+(i-1)*(n/k)):(i*(n/k))]
        
        train[i,,] <- x_split[-c((1+(i-1)*(n/k)):(i*(n/k))),]
        train_y[i,,] <-y_split[-c((1+(i-1)*(n/k)):(i*(n/k)))]
        
        
        beta_train[i,,] <- ETReg_new(train[i,,],train_y[i,,],lambda = lambda,beta.01,tau.hat[i],sigma.2[i])$coef
        
        residual[i,,] <- test_y[i,,] - cbind(rep(1,(n/k)),test[i,,])%*%beta_train[i,,]
      } 
      residual_tog <- c(residual[1:k,,])
      residual_sorted <- sort(residual_tog^2)
      residual_final <- residual_sorted[1:floor(h*(n+1))]
      # trimming the large residuals
      res[j,l] <- sqrt(mean(residual_final))
     }
    
  }
     res_avg <- colMeans(res)
    
    l <- which.min(res_avg)
    # optimaal lambda minimizing RTMSPE
    beta_solution <- beta_final[l,]
    # choosing beta with corresponding lambda
    fit=NULL
    fit$coef <- beta_solution
    # Final estimate of beta by cross-validation
    fit$lambda_opt <- l
    # optimal choice of lambda
    return(fit)
  }
  
### Performance measures (FPR, FNR, RMSPE, AUC)

ETReg_PM_sol_path = function(x_new,y_new,true_beta,beta_solution_path,beta_solution){
   # function to compute FPR, FNR,RMSPE,AUC and stepwise predictor variable 
   # selection for each simulation
   #          Args:
   #              x_new: generating a new set of design matrix 
   #              y_new: generating a new set of response using true beta
   #          true_beta: the true value of beta
   # beta_solution_path: a set of solutions of beta obtained through ETReg_sol_path
   #      beta_solution: final estimate of beta obtained thorugh ETReg_solution
 
  length <- dim(beta_solution_path)[1]
  TP_final <- rep(0,length)
  FP_final <- rep(0,length)
  p <- dim(x_new)[2]
  for(i in 1:length){
    tpr <- matrix(0,(p+1),1)
    fpr <- matrix(0,(p+1),1)
    nzero <- which(true_beta!=0)
    # non-zero elements of beta
    zero <- which(true_beta!=0)
    # zero elements of beta
    for(k1 in 1:length(nzero)){
      tpr[k1] <- as.numeric(beta_final[i,nzero[k1]]!=0) 
      # true positive numbers of predictors 
    }
    for(k2 in 1:length(zero)){
      fpr[k2] <- as.numeric(beta_final[i,zero[k2]]!=0) 
      # false positive number of predictors
    }
    TP_final[i] <- sum(tpr)/length(nzero) 
    # true positive rate
    FP_final[i] <- sum(fpr)/length(zero)
    # false positive rate
  }
  
  combind <-cbind(TP_final,FP_final)
  
  FPR <- unique(combind[,2])
  FPR <- sort(FPR) 
  TPR_initial<- rep(0,length(FPR))
  TPR <- rep(0,21)
  TPR_initial[1] <- max(unique[,1][which(unique[,2]==FPR[1])])
  for(ind in 2:length(FPR)){
    TPR_initial[ind] <- max(max(unique[,1][which(unique[,2]==FPR[ind])]),TPR_initial[ind-1])
    
  }
  
  
  for(k in 1:21)
  {
    for(g in 1:20){
    if(length(which(FPR==(k-g)*5/100)) != 0){
      TPR[k] <- TPR_initial[which(FPR==(k-g)*5/100)]
      break
    }
    }
  }
  a <- 1:21
  ROC <- cbind(TPR,(a-1)*5/100)
  # ROC, i.e., TPR corresponding to FPR
  id <- order(ROC[,2])
  
  AUC <- sum(diff(ROC[,2][id])*rollmean(ROC[,1][id],2))
  # AUC, i.e, area under the ROC curve
  
  index <- which(as.numeric(duplicated(combind))==0)-1
  index <- c(index[-1])
  
  beta_stepwise_selection <- beta_solution_path[index,] 
  # stepwise selection of predictor variables in beta
  fit <- NULL
  fit$beta_stepwise_selection <- beta_stepwise_selection
  fit$AUC <- AUC
  fit$ROC <- ROC
  # Computing FPR, FNR, RMSPE for the final beta solution
  n <- dim(x_new)
  x1 <- cbind(rep(1,n),x_new)
  RMSPE <- sqrt(mean((y-x1%*%beta_solution)^2))
  fit$RMSPE <- RMSPE
  fnr_sol <- matrix(0,(p+1),1)
  fpr_sol <- matrix(0,(p+1),1)
 
  for(k1 in length(nzero)){
    fnr_sol[k1] <- as.numeric(beta_solution[nzero[k1]]==0) 
  }
  for(k2 in length(zero)){
    fpr[k2] <- as.numeric(beta_solution[zero[k2]]!=0) 
  }
  FNR_sol <- sum(fnr_sol)/length(nzero) 
  # FNR corresponding to beta_solution
  FPR_sol <- sum(fpr)/length(zero)
  # FPR corresponding to beta_solution
  fit$FNR <- FNR_sol
  
  fit$FPR <- FPR_sol
}