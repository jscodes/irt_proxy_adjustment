###########################################################################
###                     IRT Equating within Strata                      ###
###########################################################################
# This file contains code to run the IRTeq-wStrata proxy adjustment method
# The IRT model is fit using a 2-parameter logistic IRT model

## Load libraries
library(dplyr)
library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)


## Save Stan model
#irtnoreg <- stan_model(file = "Stan Models/IRT_2PL_NoReg.stan")
#save("irtnoreg", file = "Stan Models/irtnoreg.RData")

### Load Stan model
load("Stan Models/irtnoreg.RData")


### Load example dataset
data=read.csv(file="example.csv")
#example dataset has n=500 and k=15


### Function to compute proxy adjusted values based on IRTeq-wStrata
irteq = function(k, subdata, m=50){
  
  X=subdata[,c(paste("X", 1:k, sep=""), "sum.X")]
  Y=subdata[,c(paste("Y", 1:k, sep=""), "sum.Y")]
  W=subdata$W
  Z=subdata$Z
  
  #############################################  
  ###### Run Self-Report Only Model
  #self-report patient items only data
  X.pat = X[W==0,] #for patients only
  X.pat.items = X.pat[, -(k+1)] #for all 
  X.pat.items.long <- data.frame(X.pat.items)%>%pivot_longer(cols=1:k)
  X.pat.items.long=X.pat.items.long$value
  
  #create indicators and run self-report only model
  N=length(W)-sum(W)
  P=k
  Obs <- N*P
  ii <- rep(1:N, each = P)
  kk <- rep(1:P, times = N)
  fit.pat <- sampling(irtnoreg, 
                      data = list(Obs = Obs, 
                                  N = N, 
                                  P=P, 
                                  ii=ii, 
                                  kk=kk, 
                                  x=X.pat.items.long),
                      chains=1)
  theta.pat=summary(fit.pat, pars=('theta'))$summary[,1]
  
  
  #############################################
  ###### Run Proxy-Report Only Model  
  #proxy-report patient items only data
  Y.pro = Y[W==1,] #for proxy only
  Y.pro.items = Y.pro[, -(k+1)] #for all 
  Y.pro.items.long <- data.frame(Y.pro.items)%>%pivot_longer(cols=1:k)
  Y.pro.items.long=Y.pro.items.long$value  
  
  #create indicators and run proxy only model
  N_proxy=sum(W)
  P=k
  Obs <- N_proxy*P
  ii <- rep(1:N_proxy, each = P)
  kk <- rep(1:P, times = N_proxy)
  fit.pro <- sampling(irtnoreg, 
                      data = list(Obs = Obs, 
                                  N = N_proxy, 
                                  P=P, 
                                  ii=ii, 
                                  kk=kk, 
                                  x=Y.pro.items.long),
                      chains=1)
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
  Rhat.proxy=matrix(NA,nrow=N_proxy, ncol=m)
  
  for (j in 1:m) {
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
      Rhat.proxy[v,j]= X.pat[sample(which(distance %in% sort(distance)[1:3]), 1), (k+1)]
    }
    
  }
  return(Rhat.proxy)
}



#function to run the IRT equating across all classes based on propensity score
irteq.strata.MI = function( k, data, n.class, m=50){
  
  X=data[,c(paste("X", 1:k, sep=""), "sum.X")]
  Y=data[,c(paste("Y", 1:k, sep=""), "sum.Y")]
  W=data$W
  Z=data$Z
  
  #Compute combined outcome, R
  R <- ifelse(W==1, Y$sum.Y, X$sum.X)
  
  #compute propensity score
  ps.fit=glm( W ~ Z, family="binomial")
  ps.score=predict(ps.fit, type="response")
  data$ps=ps.score
  
  #create propensity score strata as dictated by n.class
  class=cut(ps.score,
            breaks = quantile(ps.score, prob = (0:n.class)/n.class)-c(0.01 ,rep(0,n.class)), 
            labels=FALSE)
  
  #need to make sure there is at least one proxy and patient in each stratum 
  #if not, then reduce the number of classes until there are at least one in each
  n.pro<-tapply(W, class, sum)
  n.pat<-tapply((1 - W), class, sum)
  
  n.class.true=n.class  
  while((length(which(n.pro < 3)) > 0 | length(which(n.pat < 3)) > 0)) {
    
    # reduce number of class
    n.class.true=n.class.true-1
    
    # recompute the class strata with reduced number of classes
    if(n.class.true>1) {
      class = cut(ps.score,
                  breaks = quantile(ps.score, prob = (0:n.class.true)/n.class.true)-c(0.01 ,rep(0,n.class.true)), 
                  labels=FALSE)
    } else {
      class = rep(1, length(ps.score))
    }
    n.pro<-tapply(W, class, sum);
    n.pat<-tapply((1 - W), class, sum);  
  }
  
  #Compute combined adjusted outcome, R.hat, by class
  R.hat<- matrix(NA, nrow=length(R), ncol=m)  
  
  # adjust the proxy (source) values
  for(c in 1:n.class.true) {
    R.hat[class==c & W==1,]= irteq(k, data[class==c,]) 
  }
  
  # the patient (destination) values remain as observed 
  R.hat[W==0,] = R[W==0]
  
  save.Rhat=as.list(as.data.frame(R.hat))

  return(save.Rhat)  
}

### Run on example dataset
# R.hat.irteqstr= irteq.strata.MI(k=15, data, m=50, n.class=5)

















    

      



