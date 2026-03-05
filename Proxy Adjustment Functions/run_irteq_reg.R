###########################################################################
###             IRT Equating with Regression Adjustment                 ###
###########################################################################
# This file contains code to run the IRTeq-Reg proxy adjustment method
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



### Function to compute proxy adjusted values based on IRTeq-Reg
irteq.reg.MI = function(k, data, m=50){
  
  X=data[,c(paste("X", 1:k, sep=""), "sum.X")]
  Y=data[,c(paste("Y", 1:k, sep=""), "sum.Y")]
  W=data$W
  Z=data$Z
  
  #Compute combined outcome, R
  R <- ifelse(W==1, Y$sum.Y, X$sum.X)
  
  #############################################  
  ###### Run Self-Report Only Model
  #self-report patient items only data
  X.pat = X[W==0,] #for patients only
  X.pat.items = X.pat[, -(k+1)] #for all 
  X.pat.items.long <- data.frame(X.pat.items)%>%pivot_longer(cols=1:k)
  X.pat.items.long=X.pat.items.long$value
  
  #standardize covars
  Z.pat=Z[W==0]
  Z.pat.std=(Z.pat-mean(Z.pat))/sd(Z.pat)
  
  #create indicators and run self-report only model
  N=length(W)-sum(W)
  P=k
  Obs <- N*P
  ii <- rep(1:N, each = P)
  kk <- rep(1:P, times = N)
  fit.pat <- sampling(irtreg1, 
                      data = list(Obs = Obs, 
                                  N = N, 
                                  P=P, 
                                  ii=ii, 
                                  kk=kk, 
                                  z=Z.pat.std, 
                                  x=X.pat.items.long),
                      chains=1, iter=1000)
  theta.pat=summary(fit.pat, pars=('theta'))$summary[,1]
  
  
  #############################################
  ###### Run Proxy-Report Only Model  
  #proxy-report patient items only data
  Y.pro = Y[W==1,] #for proxy only
  Y.pro.items = Y.pro[, -(k+1)] #for all 
  Y.pro.items.long <- data.frame(Y.pro.items)%>%pivot_longer(cols=1:k)
  Y.pro.items.long=Y.pro.items.long$value  
  
  #standardize covars
  Z.proxy=Z[W==1]
  Z.proxy.std=(Z.proxy-mean(Z.proxy))/sd(Z.proxy)  
  
  #create indicators and run proxy only model
  N_proxy=sum(W)
  P=k
  Obs <- N_proxy*P
  ii <- rep(1:N_proxy, each = P)
  kk <- rep(1:P, times = N_proxy)
  fit.pro <- sampling(irtreg1, 
                      data = list(Obs = Obs, 
                                  N = N_proxy, 
                                  P=P, 
                                  ii=ii, 
                                  kk=kk, 
                                  z=Z.proxy.std, 
                                  x=Y.pro.items.long),
                        chains=1, iter=1000)
  theta.pro.means=summary(fit.pro, pars=('theta'))$summary[,c(1)]  
  theta.pro.sd=summary(fit.pro, pars=('theta'))$summary[,c(3)]
  
  #############################################
  #####   Compute CDFs for thetas using KDE
  
  #self-report cdf
  theta.pat.density=density(theta.pat) 
  cdf.pat=data.frame(cbind(theta.vals=theta.pat.density$x, 
                           cumprob=cumsum(theta.pat.density$y)/sum(theta.pat.density$y)))   
  
  
  #############################################
  ##### Equate proxy outcome values to self-report cdf
  # create vector to save imputed values
  
  #Save multiply imputed datasets to output
  save.Rhat=list()
  
  #v=1
  for (j in 1:m) {
  
  Rhat.proxy=rep(NA,N_proxy)
  
    #v=1
    for (v in 1:N_proxy){
      
      #resample proxies
      theta.pro=mapply(rnorm, n=rep(1, length(theta.pro.means)), mean=theta.pro.means, sd=theta.pro.sd)
      
      #compute the proxy-report cdf
      theta.pro.density=density(theta.pro) 
      cdf.pro=data.frame(cbind(theta.vals=theta.pro.density$x, 
                               cumprob=cumsum(theta.pro.density$y)/sum(theta.pro.density$y))) 
      
      #find the percentile for each proxy-report, then find the self-report with the closest
      #percentile, and the used the observed value of the self-report as the imputed value for
      #the proxy-reports outcome value
      
      #percentile for proxy
      proxy.p =cdf.pro[which.min(round(abs(cdf.pro[,1]-theta.pro[v]),4)),2]    
      
      #find corresponding self-report cdf ability of percentile (proxy-percentile)
      equated.theta=cdf.pat[which.min(round(abs(cdf.pat[,2]-proxy.p),4)),1]
      
      #find observed outcome value of closest equated ability self-report ability and impute as rhat
      distance=round(abs(theta.pat-equated.theta),4)
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
# R.hat.irteqreg= irteq.reg.MI(k=15, data, m=50)





















    

      



