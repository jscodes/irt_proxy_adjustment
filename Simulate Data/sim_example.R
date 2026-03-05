
###########################################################################
###                CODE TO SIMULATE PROXY-REPORTED DATA                 ###
###########################################################################

## Load libraries
library(dplyr)
library(tidyverse)
library(MASS)

########################################################################
###                 Necessary data creation functions                ###
########################################################################


####################################
### Function to take inverse logit(x)
expit<-function(x) exp(x)/(1+exp(x))


####################################
### Function to compute alpha for shifting the intercept to have 
#proportion 'missing' (ie, proxies) as mean
fun.alpha<-function(covar, b, mean) {
  fun<-function(alpha) {
    {mean(expit(alpha + b*covar)) - mean}}
  return(uniroot(fun, c(-100,100)))
}


####################################
### Functions that find the values of psi for a specific correlation
find.delta = function(cor) {
  fun.corXZ=function(delta){ delta/sqrt(1+delta^2)-cor}
  return(uniroot(fun.corXZ, c(-100,100))$root)
}

find.psi.x = function(cor, delta, psi.z) {
  fun.corXY.x=function(psi.x){ 
    (delta*(delta+psi.z+delta*psi.x) + (1+psi.x)  ) / 
      (sqrt(delta^2+1)*sqrt((delta+psi.z+delta*psi.x)^2 +(1+psi.x)^2 +1)  ) -cor}
  return(uniroot(fun.corXY.x, c(-100,100))$root)
}

find.psi.z = function(cor, delta, psi.x) {
  fun.corXY.z=function(psi.z){ 
    (delta*(delta+psi.z+delta*psi.x) + (1+psi.x)  ) / 
      (sqrt(delta^2+1)*sqrt((delta+psi.z+delta*psi.x)^2 +(1+psi.x)^2 +1)  ) -cor}
  return(uniroot(fun.corXY.z, c(-100,0))$root)
}


######################################################################
###         SIMULATE SELF-REPORT, X, and PROXY-REPORT, Y           ###
######################################################################
# Simulate from a multivariate norm and then discretize based on 
# randomly drawn cutoff values from uniform dist


### Function to simulate one data set from MVN model
sim.MVN.xy = function(n, k, z, delta, sigma, psi.z, psi.x, offset=0) {
  
  # normalized errors
  corrmatrix=matrix(sigma, nrow=k, ncol=k)
  diag(corrmatrix)=1
  norm_err=mvrnorm(n, rep(0,k), Sigma = corrmatrix) 
  
  # normalized responses
  x_ik_norm = delta*z + norm_err
  
  # proxy responses: X plus measurement error based on Z or X
  y_ik_norm=x_ik_norm  + psi.z*z + psi.x*x_ik_norm + rnorm(n,1)
  
  #generate random values for proportion of yes responses
  #props=(1/(k+1))*seq(1,k)
  props=runif(k, .3, .7)
  x_ik=matrix(NA, nrow=n,ncol=k)
  y_ik=matrix(NA, nrow=n,ncol=k)
  cutoff=ifelse(props+offset<0.01,0.01,props+offset)
  cutoff=ifelse(cutoff>0.99,0.99,cutoff)    
  for (i in 1:k)
  {
    x_ik[,i] = 1*( scale(x_ik_norm[,i]) >= qnorm(props[i]) ) 
    y_ik[,i] = 1*( scale(y_ik_norm[,i]) >= qnorm(cutoff[i]) ) 
  }
  
  #evaluate X_ik
  sum.x_i = rowSums(x_ik)
  sum.y_i = rowSums(y_ik)  
  
  #X data
  Xdata=data.frame(cbind(x_ik, sum.x_i))
  colnames(Xdata)=c(paste("X", 1:k, sep=""), "sum.X")
  
  #Y data
  Ydata=data.frame(cbind(y_ik, sum.y_i))
  colnames(Ydata)=c(paste("Y", 1:k, sep=""), "sum.Y")
  
  return(list(X=Xdata, Y=Ydata))
}



######################################################################
###                     SIMULATE PROXY-INDICATOR, W                ###
######################################################################


####################################
### Function to compute prob of being proxy and proxy indicator, W
sim.proxies = function(miss.covar, gamma, mu ) {
  alpha=fun.alpha(miss.covar, gamma, mean=mu)$root 
  p.W= expit(alpha+gamma*miss.covar) #mean for predicted probabilities is 0.5 with alpha=0
  W=rbinom(n=length(miss.covar), size=1, prob=p.W )
  return(W)
}





########################################################################
###                     EXAMPLE SIMULATED DATA                       ###
########################################################################


### Function for simulating a full dataset with W, Z, X, Y, R
sim.dataset = function(
    n, #number of subjects
    mu, #proportion of proxies
    k, #number of items
    miss.covar, #missingness mechanism covariate
    gamma, #strength of missingness mechanism
    offset, #omega
    delta.cor,
    psi.z.cor,
    psi.x.cor)  {
  
  d=find.delta(delta.cor)
  
  #Simulate covariate, Z
  Z=scale(rnorm(n))[,1]
  
  if (psi.z.cor != 0){
    XYdata=sim.MVN.xy(n=n, 
                       k=k, 
                       z=Z, 
                       delta=d, 
                       sigma=.3, 
                       psi.z=find.psi.z(psi.z.cor, delta=d, psi.x=0), 
                       psi.x=0,  
                       offset=offset) 
  } else if (psi.x.cor != 0) {
    XYdata=sim.MVN.xy(n=n, 
                       k=k, 
                       z=Z, 
                       delta=d, 
                       sigma=.3, 
                       psi.z=0, 
                       psi.x=find.psi.x(psi.x.cor, delta=d, psi.z=0),  
                       offset=offset) 
  } else {
    XYdata=sim.MVN.xy(n=n, 
                       k=k, 
                       z=Z, 
                       delta=d, 
                       sigma=.3, 
                       psi.z=0, 
                       psi.x=0,  
                       offset=offset) 
  }
  
  
  X=XYdata$X
  Y=XYdata$Y
  data=cbind(Z,X,Y)
  
  #Simulate proxy indicator, W
  data$W=sim.proxies(miss.covar=data[,miss.covar], gamma=gamma, mu=mu)  
  
  #Operationalize observed data R
  R=X
  R[W==1,] = Y[W==1,]
  colnames(R)=c(paste("R", 1:k, sep=""), "sum.R")
  
  
  #put all data together in frame
  data.all=data.frame(cbind(data,R))   
  
  #save simulated dataset
  return(data.all)
  
}
  
  
  

## parameter values use for example dataset
n=c(500)
mu=c(0.10)
k=c(15)
miss.covar=c("Z") #c("sum.X")
gamma=c(0.5)
offset=c(0)
delta.cor=c(0.3)
psi.z.cor=c(0)
psi.x.cor=c(0.3)

set.seed(23)
example.data= sim.dataset(n, 
                          mu, 
                          k, 
                          miss.covar, 
                          gamma, 
                          offset, 
                          delta.cor,
                          psi.z.cor,
                          psi.x.cor)

#write.csv(example.data, file="example.csv")




