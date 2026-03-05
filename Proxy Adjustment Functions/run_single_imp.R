###########################################################################
###                      SINGLE IMPUTATION METHODS                      ###
###########################################################################
# This file contains code for the single imputation methods:
# Substitution
# Regression Adjustment
# Regression Prediction


## Load libraries
library(dplyr)
library(tidyverse)

### Load example dataset
data=read.csv(file="example.csv")
#example dataset has n=500 and k=15


###########################################################
###                SUBSTITUTION FUNCTION                ###
# no adjustment or imputation just use substituted data and 

substitution = function(data){
  
  #Pull sum.R
  R.hat = data$sum.R
  
  return(R.hat)  
}
## Run on example dataset
# R.hat.sub=substitution(data)



###########################################################
###                REGRESSION ADJUSTMENT                ###
# regress on combined outcome var, R scores with covariates 
# and an indicator for proxy response, W and use fitted values 
# for proxies assuming self-report responses


reg.adjust = function(k, #number of items
                      data
                      ){
  
  #fit regression model with covariate and W as a predictor
  reg.fit=lm(sum.R ~ Z+W, data=data)

  #compute design matrix for patients with proxies
  Xmat.pred = data.frame(cbind(model.matrix(~Z, data=data)))
  Xmat.pred$W=0 #predict as self-report
  Xmat.pred[,c("pred")] =predict(reg.fit, Xmat.pred)
  
  #Compute combined adjusted outcome, R.hat, round to the nearest whole number
  R.hat <- round(ifelse(data$W==1, Xmat.pred$pred, data$sum.R), 0)
  
  #put min and max values for predicted value to avoid extrapolation
  R.hat <- ifelse(R.hat<0, 0, R.hat)
  R.hat <- ifelse(R.hat>p, p, R.hat)
  
    
  return(R.hat) 
}

## Run on example dataset
# R.hat.regadj=reg.adjust(15, data)




###########################################################
###                REGRESSION PREDICTION                ###
# regress on patient scores only using covariate, Z and then 
# predict proxy-reported data using linear model estimates

reg.predict = function(k, #number of items
                       data
                       ){
  
  #fit regression model with covariate on only patient responses
  reg.fit=lm(sum.X ~ Z, data=data[data$W==0,])
  
  #summary(reg.fit)
  Xmat.pred = data.frame(cbind(model.matrix(~Z, data=data)))
  Xmat.pred[,c("pred")] =predict(reg.fit, Xmat.pred)
  
  #Compute combined adjusted outcome, R.hat,round to the nearest whole number
  R.hat <- round(ifelse(data$W==1, Xmat.pred$pred, data$sum.R), 0)
  
  #put min and max values for predicted value to avoid extrapolation
  R.hat <- ifelse(R.hat<0, 0, R.hat)
  R.hat <- ifelse(R.hat>p, p, R.hat)
  

  
  return(R.hat) 
}

## Run on example dataset
# R.hat.regpred=reg.predict(15, data)







    

      



