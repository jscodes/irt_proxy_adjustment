
###########################################################################
###                IRT-Latent Regression (all patients)                 ###
###########################################################################
# This file contains code to run the IRT-Reg2 proxy adjustment method
# The IRT model is fit using a 2-parameter logistic IRT model

## Load libraries
library(dplyr)
library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


## Save Stan model
#irtreg2 <- stan_model(file = "Stan Models/IRT_2PL_Reg2.stan")
#save("irtreg2", file = "Stan Models/irtreg2.RData")

### Load Stan model
load("Stan Models/irtreg2.RData")


### Load example dataset
data=read.csv(file="example.csv")
#example dataset has n=500 and k=15


### Function to compute proxy adjusted values based on IRT-Reg1
irt.reg2.MI = function(k, data, m=50){
  
  X=data[,c(paste("X", 1:k, sep=""), "sum.X")]
  Y=data[,c(paste("Y", 1:k, sep=""), "sum.Y")]
  W=data$W
  Z=data$Z
  
  #Compute combined outcome, R
  R <- ifelse(W==1, Y$sum.Y, X$sum.X) 
  
  # X values (true/self-report) long form by item
  X.pat = X[W==0,] 
  X.items = X[, -(k+1)] #for all 
  X.items.long <- data.frame(X.items)%>%pivot_longer(cols=1:k)
  X.items.long=X.items.long$value
  
  # Y values (proxy report) long form by item
  Y.items = Y[, -(k+1)] #for all 
  Y.items.long <- data.frame(Y.items)%>%pivot_longer(cols=1:k)
  Y.items.long=Y.items.long$value  
  
  # W value long form
  W.rep=matrix(W, nrow =length(W), ncol=k,byrow=F)
  W.long=data.frame(W.rep)%>%pivot_longer(cols=1:(k))
  W.long=W.long$value
  
  # R value long form (combined outcome)
  R.items.long <-ifelse(W.long==1, Y.items.long, X.items.long)
  
  
  #standardize covars
  Z.std=(Z-mean(Z))/sd(Z)
  
  #create indicators
  N=length(W)
  P=k
  Obs <- N*P
  ii <- rep(1:N, each = P)
  kk <- rep(1:P, times = N)
  
  
  fit <- sampling(irtreg2, 
                  data = list(Obs = Obs, 
                              N = N, 
                              P=P, 
                              ii=ii, 
                              kk=kk,
                              z=Z.std, 
                              w=W, 
                              x=R.items.long),
                  chains=1)

  #predicted means for patients
  params=summary(fit, pars=c("delta","rho"))$summary[c(1,2),c(1,3)] 
  delta = params[1,1]
  rho = params[2,1]
  pred.theta.pat=delta*Z.std[W==0]#+rho*W[W==0]
  
  #Save multiply imputed datasets to output
  save.Rhat=list()
  
  #v=1
  for (j in 1:m) {
    
    #Impute R values for proxies
    Rhat.proxy=rep(NA,sum(W))
    #v=1
    for (v in 1:sum(W)){
      
      #take a random draw of delta to compute proxy values
      delta=rnorm(1,params[1,1], params[1,2])
      rho=rnorm(1,params[2,1], params[2,2])
      pred.theta.pro=delta*Z.std[W==1]#+rho*W[W==1]

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
# R.hat.irtreg2=irt.reg2.MI(k=15, data, m=50)





