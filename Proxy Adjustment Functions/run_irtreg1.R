###########################################################################
###            IRT-Latent Regression (self-report only)                 ###
###########################################################################
# This file contains code to run the IRT-Reg1 proxy adjustment method
# The IRT model is fit using a 2-parameter logistic IRT model

## Load libraries
library(dplyr)
library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)

## Save Stan model
#irtreg1 <- stan_model(file = "Stan Models/IRT_2PL_Reg1.stan")
#save("irtreg1", file = "Stan Models/irtreg1.RData")

### Load Stan model
load("Stan Models/irtreg1.RData")


### Load example dataset
data=read.csv(file="example.csv")
#example dataset has n=500 and k=15


### Function to compute proxy adjusted values based on IRT-Reg1
irt.reg1.MI = function(k, data, m=50){
  
  X=data[,c(paste("X", 1:k, sep=""), "sum.X")]
  Y=data[,c(paste("Y", 1:k, sep=""), "sum.Y")]
  W=data$W
  Z=data$Z
  
  #Compute combined outcome, R
  R <- ifelse(W==1, Y$sum.Y, X$sum.X) 
  
  
  # X values (true/self-report) long form by item
  X.pat = X[W==0,]
  X.pat.items = X.pat[, -(k+1)] #for all 
  X.pat.items.long <- data.frame(X.pat.items)%>%pivot_longer(cols=1:k)
  X.pat.items.long=X.pat.items.long$value
  
  
  #standardize covars
  Z.pat=Z[W==0]
  Z.pat.std=(Z.pat-mean(Z.pat))/sd(Z.pat)
  Z.proxy=Z[W==1]
  Z.proxy.std=(Z.proxy-mean(Z.pat))/sd(Z.pat)  
  
  #create indicators
  N=length(W)-sum(W)
  N_proxy=sum(W)
  P=k
  Obs <- N*P
  ii <- rep(1:N, each = P)
  kk <- rep(1:P, times = N)
  irt2P.fit <- sampling(irtreg1, 
                        data = list(Obs = Obs, 
                                    N = N, 
                                    P=P, 
                                    ii=ii, 
                                    kk=kk,
                                    z=Z.pat.std, 
                                    x=X.pat.items.long),
                        chains=1)
  
  #predicted means for patients
  params=summary(irt2P.fit, pars=c("delta"))$summary[1,c(1,3)] 
  delta = params[1]
  pred.theta.pat=delta*Z.pat.std
  
  #Save multiply imputed datasets to output
  save.Rhat=list()
  
  #v=1
  for (j in 1:m) {
    
    #Impute R values for proxies
    Rhat.proxy=rep(NA,N_proxy)
    #v=1
    for (v in 1:N_proxy){
      #take a random draw of delta to compute proxy values
      delta=rnorm(1,params[1], params[2])
      pred.theta.pro=delta*Z.proxy.std

      #take random draws of theta for patients
      #finding the 3 closest abilities of the patient data for each proxy,
      #then randomly sample observed value from one of them if there are ties to impute proxy value 
      distance= abs(pred.theta.pat-pred.theta.pro[v]) 
      Rhat.proxy[v]= X.pat[sample(which(distance %in% sort(distance)[1:3]), 1), (k+1)]
    }  
  
  #Compute combined adjusted outcome, R.hat
  R.hat<-R  
  R.hat[which(W==1)] <- Rhat.proxy
  
  save.Rhat[[j]]=R.hat
  
}
  
  return(save.Rhat)  
}


### Run on example dataset
# R.hat.irtreg1=irt.reg1.MI(k=15, data, m=50)









    

      



