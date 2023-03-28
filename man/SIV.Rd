\name{SIV}
\alias{SIV}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  the proposed sythetic instrument variable estimator for causal effect estimation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SIV(X,Y,q = "NULL")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{    the n*p design matrix    }
  \item{Y}{    the n*1 outcome vector    }
  \item{q}{    number of latent confounder       we use onatski (2009)'s method as default method }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
 \item{beta}{the solution of $l0$ optimization}
 \item{identifiablity}{model identifiability}
 \item{Xhat}{first stage result after projection}
 \item{qhat}{estimated number of confounder}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
  Dingke Tang
 dingke.tang@mail.utoronto.ca
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
p = 100
s = 5
q = 2
n = 500
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
}
