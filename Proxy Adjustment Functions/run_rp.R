###########################################################################
###                         RANK PERMUTATION                            ###
###########################################################################
# This file contains code to run the RP proxy adjustment method


### Load libraries
library(dplyr)
library(tidyverse)


### Load example dataset
data=read.csv(file="example.csv")
#example dataset has n=500 and k=15


### Function to compute proxy adjusted values based on RP
#1) fit a logistic model predicting probability of being proxy adjusted by Z,
#2) create subgroups of n patients based on propensity of being proxy
#3) run RP within each of the subgroups
RP = function(subdata, min.score, max.score){
  
  W=subdata$W
  R=subdata$sum.R
  
  dest=R[W==0] #patients
  source=R[W==1] #proxies
  
  # Compute ranks
  dest.rank=rank(dest, ties.method = "random")
  source.rank=rank(source, ties.method ="random")
  
  # Combine the source and destination samples together, then sample ordering
  total=c(rep(1,length(dest)),rep(0,length(source)))
  total.order=sample(total,replace=F)
  
  # Extract overall ranks for destination and source 
  total.rank=c(1:length(total))
  total.dest.rank=total.rank[total.order==1]
  total.source.rank=total.rank[total.order==0]  
  
  #Combine respective source and destination ranks for total sample
  combined.rank=rep(0,length(total))
  combined.rank[total.order==1]=sort(dest.rank)
  combined.rank[total.order==0]=sort(source.rank)
  
  ## Imputing values on the Destination scale
  total.dest.score=rep(NA,length(total))
  
  for (i in 1:length(total.dest.score)){
    
    #Destination scores retain their original value
    if (total.order[i]==1){
      total.dest.score[i]=dest[dest.rank==combined.rank[i]]}
    
    #Source scores are randomly sampled from Gaussian kernel density
    else{
      #If Source rank is the lowest sampled pooled rank, set min possible value to min of Destination scores
      if (i==1){
        lower.bound=min.score
        upper.bound=min(dest)
        dest.density=density(dest,kernel=c("gaussian"),from=lower.bound,to=upper.bound)$x
        total.dest.score[i]=sample(dest.density,1)
      }
      #If Source rank is the highest sampled rank, set max possible value to max of Destination scores
      else if (i==length(total.dest.score)){
        lower.bound=max(dest)
        upper.bound=max.score
        dest.density=density(dest,kernel=c("gaussian"),from=lower.bound,to=upper.bound)$x
        total.dest.score[i]=sample(dest.density,1)
      }
      
      else{
        lower.bound=total.dest.score[i-1]
        
        upper.bound=suppressWarnings(dest[dest.rank==min(combined.rank[total.rank>i & total.order==1])])
        if (length(upper.bound)==0) { upper.bound=max.score}
        
        dest.density=density(dest,kernel=c("gaussian"),from=lower.bound,to=upper.bound)$x
        total.dest.score[i]=sample(dest.density,1)
      }
    }
  }
  
  #Rank permuted source scores
  source.vals=rep(0,length(source))
  for (i in 1:length(source.vals)){
    source.vals[i]=total.dest.score[total.order==0 & combined.rank==source.rank[i]]
  }
  source.vals
  return(source.vals)
}


RPinClass = function(k, data, n.class, min.score, max.score){
  
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
  
  #create propensity score strat as dictated by n.class
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
  R.hat<- rep(NA, length(R))  
  
  # adjust the proxy (source) values
  for(c in 1:n.class.true) {
    R.hat[class==c & W==1]= round(RP(data[class==c,], min.score, max.score),0)
  }
  
  # the patient (destination) values remain as observed 
  R.hat[W==0] = R[W==0]
  
  
  return(R.hat)  
}


### Run on example dataset
# R.hat.rp=RPinClass(k=15, data=data, n.class=5, min.score=0, max.score=15)























    

      



