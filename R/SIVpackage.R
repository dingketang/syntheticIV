if(!require("MASS")) install.packages("MASS")
if(!require("glmnet")) install.packages("glmnet")
if(!require("abess")) install.packages("abess")
if(!require("caret")) install.packages("caret")
if(!require("leaps")) install.packages("leaps")

library(MASS)
library(glmnet)
library(abess)
library(caret)
library(leaps)


Tuning_q<-function(Data)
{
  p<-ncol(Data)
  n<-nrow(Data)
  rmax<-min(10,(p-1)/2)
  X=scale(Data,center = TRUE,scale = FALSE)
  eigenv<-eigen(cov(X))$values
  teststat <- rep(1,rmax)
  for (i in 1:rmax) {
    teststat[i] = (eigenv[i] - eigenv[i+1])/(eigenv[i+1] - eigenv[i+2])
  }
  which.max(teststat)
}


get_model_formula <- function(id, object, outcome){
  models <- summary(object)$which[id,-1]
  predictors <- names(which(models == TRUE))
  predictors <- paste(predictors, collapse = "+")
  as.formula(paste0(outcome, "~", predictors))
}

get_cv_error <- function(model.formula, data){
  set.seed(1)
  train.control <- trainControl(method = "cv", number = 10)
  cv <- train(model.formula, data = data, method = "lm",
              trControl = train.control)
  cv$results$RMSE
}



##########################################
# Description: exact solution of (8)     #
# Input:                                 #
#  X: the n*p design matrix              #
#  Y: the n*1 outcome vector             #
#  q: number of latent confounder        #
# Output:                                #
# beta: the solution of l_0 optimization #
##########################################

CV_Best_SS<- function(X,Y, q){
  n = nrow(X)
  p = ncol(X)
  dataset <- data.frame(Y = Y , X  = X)
  names(dataset) = c("Y",paste0("X",1:p))

  models <- regsubsets(Y~.,data = dataset, nvmax = p-q)
  cv.errors <- rep(0,p-q)
  for (i in 1:(p-q)) {
    cv.errors[i]=get_cv_error(get_model_formula(i, models, "Y"),dataset)
  }
  betaresult <- rep(0,p)
  names(betaresult) <-paste0("X",1:p)
  betaresult[names(coef(models,which.min(cv.errors))[-1])]=
    coef(models,which.min(cv.errors))[-1]
  betaresult
}

##############################################
# Description: the propsed estimator         #
# Input:                                     #
#  X: the n*p design matrix                  #
#  Y: the n*1 outcome vector                 #
#  q: number of latent confounder            #
#  we use onatski(2009)'s method as default  #
#  method to find number of latent confounder#
# Output:                                    #
# beta: the solution of l_0 optimization     #
# identifiablity: model identifiability      #
# Xhat : first stage result after projection #
# qhat : estimated number of confounder      #
##############################################

SIV<- function(X,Y,q = NULL){
  if(is.null(q))
  {
    q = as.numeric(Tuning_q(X))
  }
  X = scale(X,center = TRUE,scale = FALSE)
  p = ncol(X)
  n = nrow(X)
  if(n<=p){# highd method
    factor_estimate <-prcomp(X)
    if(q>1){
      load <- factor_estimate$rotation[,1:q]%*%
        diag(factor_estimate$sdev[1:q])
    }
    else
      load <- matrix(factor_estimate$rotation[,1]*
                       (factor_estimate$sdev[1]))
  }

  else{# low-d case
    fac  = factanal(X,factors = as.numeric(q))
    load = diag(sqrt(diag(var(X))))%*% fac$loadings
  }

  B_orthogonal <-  Null(load)

  X_SIV = X%*%B_orthogonal
  X_hat = fitted(lm(X~X_SIV))
  timeb = Sys.time()


  abess_fit <- abess(X_hat, Y,tune.type = "cv",nfold = 10)

  if(p>30){
    result = coef(abess_fit,
                  support.size =
                    abess_fit$support.size[which.min(abess_fit$tune.value[which(abess_fit$support.size <= (p-q))])])[-1]
  }
  else{
    result = CV_Best_SS(X_hat,Y,q)
  }
  s = sum(result!=0)
  result = list(betahat = result)
  if(s+1+q <= p) result$identifiablity = TRUE
  else           result$identifiablity = FALSE
  result$qhat = q
  result$Xhat = X_hat
  result
}
