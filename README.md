# syntheticIV
the code and package for paper ''The synthetic instrument: From sparse association to sparse causation''


## Download
The package can be downloaded in R with
```
devtools::install_github("dingketang/syntheticIV")
```

## Example
Here is a basic example illustrating the usage of the package.


### Geneate data

```
p = 100 #number of variables 
s = 5   #number of causal variables  
q = 2   #number of latent confounder
n = 500 #number of samples 
L = matrix(runif(p*q,max = 1,min = -1),nrow = p)# loading matrix
gamma = matrix(runif(q,max = 1,min = -1))       # gamma matrix
beta  = matrix(rep(0,p))
beta[1:s] = 1
U = matrix(rnorm(n*q),nrow = n)
X = U%*%t(L) + matrix(rnorm(n*p,sd = 2),nrow = n)
Y = X%*%beta + U%*%gamma + rnorm(n)*1
```

### Call function
```
SIVfit <- syntheticIV::SIV(X,Y)
SIVfit$beta
```
