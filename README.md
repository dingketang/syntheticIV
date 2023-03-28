# syntheticIV
the code and package for ''The synthetic instrument: From sparse association to sparse causation''


example:

devtools::install_github("dingketang/syntheticIV")
library("syntheticIV")

p = 100 #number of variables 
s = 5   #number of causal variables  
q = 2   #number of latent confounder
n = 500 #number of samples 
set.seed(20221224)
L = matrix(runif(p*q,max = 1,min = -1),nrow = p)# loading matrix
gamma = matrix(runif(q,max = 1,min = -1))       # gamma matrix
beta  = matrix(rep(0,p))
beta[1:s] = 1
U = matrix(rnorm(n*q),nrow = n)
X = U\%*\%t(L) + matrix(rnorm(n*p,sd = 2),nrow = n)
Y = X\%*\%beta + U\%*\%gamma + rnorm(n)*1
SIVfit <- SIV(X,Y)
SIVfit$beta
