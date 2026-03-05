###########################################################################
###                     PREDICTIVE MEAN MATCHING                        ###
###########################################################################
# This file contains code to run the PMM proxy adjustment method

## Load libraries
library(dplyr)
library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


## Save Stan model
#linear <- stan_model(file = "Stan Models/SimpleLinearReg.stan")
#save("linear", file = "Stan Models/linear.RData")

### Load Stan model
load("Stan Models/linear.RData")


### Load example dataset
data=read.csv(file="example.csv")
#example dataset has n=500 and k=15


### Function to compute proxy adjusted values based on PMM
pmmMI = function(k, data, m=50){
  
  X=data[,c(paste("X", 1:k, sep=""), "sum.X")]
  Y=data[,c(paste("Y", 1:k, sep=""), "sum.Y")]
  W=data$W
  Z=data$Z
  
  #Patient only data
  X.pat = X[W==0,] 
  
  #Compute combined outcome, R
  R <- ifelse(W==1, Y$sum.Y, X$sum.X)
  
  #Run linear model on self-reported patients only
  N=length(W)-sum(W) #number of self-reports
  N_proxy=length(W)-N #number of proxies
  
  #Standardize covars
  Z.pat=Z[W==0]
  Z.pat.std=(Z.pat-mean(Z.pat))/sd(Z.pat)
  Z.proxy=Z[W==1]
  Z.proxy.std=(Z.proxy-mean(Z.pat))/sd(Z.pat)

  #Use Stan to fit linear model
  pmm.fit <- sampling(linear, 
                      data = list(N = N, 
                                  z=Z.pat.std, 
                                  x=X$sum.X[W==0]),
                      chains=1,
                      iter = 2000)
  
  #Alternately, can use OLS linear model for MLE estimates
  #summary(lm(X.pat$sum.X~Z.pat.std))   
  
  #Compute predicted means for patients
  params=summary(pmm.fit)$summary[1:3,c(1,3)] #mean and sd for normal dist
  alpha = params[1,1]
  beta = params[2,1]
  pred.x.pat=alpha+beta*Z.pat.std

  #Save multiply imputed datasets to output
  save.Rhat=list()
  
  #j=1
  for (j in 1:m) {
    #Impute R values for proxies
    Rhat.proxy=rep(NA,sum(W))
    #v=1
    for (v in 1:N_proxy){
      #resample proxies for each match
      alpha=rnorm(1,params[1,1], params[1,2])
      beta=rnorm(1,params[2,1], params[2,2])
      pred.x.pro=alpha+beta*Z.proxy.std

      #finding the 3 closest ps of the patient data for each proxy,
      #then randomly sample observed value from one of them to impute proxy value 
      distance= abs(pred.x.pat-pred.x.pro[v]) 
      Rhat.proxy[v]= X.pat[sample(which(distance %in% sort(distance)[1:3]), 1), (k+1)]
    }   
  
    #Compute metrics
    
    #Compute combined adjusted outcome, R.hat
    R.hat<-R  
    R.hat[which(W==1)] <- Rhat.proxy
    
    save.Rhat[[j]]=R.hat
    
  }
  
  return(save.Rhat)   
}




### Run on example dataset
# R.hat.pmm=pmmMI(k=15, data, m=50)
















    

      



