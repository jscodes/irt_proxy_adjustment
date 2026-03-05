###########################################################################
###                PROPENSITY SCORE MULTIPLE IMPUTATION                 ###
###########################################################################
# This file contains code to run the PSMI proxy adjustment method

## Load libraries
library(dplyr)
library(tidyverse)

### Load example dataset
data=read.csv(file="example.csv")
#example dataset has n=500 and k=15


#Function to compute proxy adjusted values based on PSMI
pscoreMI = function(k, data, m=50, n.class=5){
  
  X=data[,c(paste("X", 1:k, sep=""), "sum.X")]
  Y=data[,c(paste("Y", 1:k, sep=""), "sum.Y")]
  W=data$W
  Z=data$Z
  
  
  #Compute combined outcome, R
  R <- ifelse(W==1, Y$sum.Y, X$sum.X)
  
  #Compute propensity score
  ps.fit=glm( W ~ Z, family="binomial")
  ps.score=predict(ps.fit, type="response")
  
  #Create propensity score strata as dictated by n.class
  class=cut(ps.score,
            breaks = quantile(ps.score, prob = (0:n.class)/n.class)-c(0.01 ,rep(0,n.class)), 
            labels=FALSE)
  
  #Make sure there is at least 3 proxy and patient in each stratum 
  #if not, then reduce the number of classes until there are at least 3 in each
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
  
  #Save multiply imputed datasets to output
  save.Rhat=list()
  
  #j=1
  for (j in 1:m) {
    
    #Compute combined adjusted outcome, R.hat, by class
    R.hat<- rep(NA, length(R))  
    
    # adjust the proxy (source) values
    #c=1
    for(c in 1:n.class.true) {
      n.pat=length(R[class==c & W==0])
      n.proxy=length(R[class==c & W==1])
      
      donors=sample(R[class==c & W==0], n.pat, replace=T)
      proxy.imputations = sample(donors, n.proxy, replace=T)
      R.hat[class==c & W==1]= proxy.imputations
    }    
    
    
    # Self-report patient values remain as observed 
    R.hat[W==0] = R[W==0]
    
    # Save adjusted values
    save.Rhat[[j]]=R.hat
  
  }
  
  return(save.Rhat)  
}



## Run on example dataset
# R.hat.psm=pscoreMI(k=15, data, m=50, n.class=5)



















    

      



