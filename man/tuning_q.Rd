\name{Tuning_q}
\alias{Tuning_q}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
the default method in SIV package to estimate number of latent confounder
}
\description{
   onatski(2009)'s method estimate number of latent confounder
}
\usage{
Tuning_q(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{the n*p design matrix X  with n samples and p variable}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
   return the number of latent confounder
}
\references{
Testing Hypotheses About the Number of Factors in Large Factor Models
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
q = 2
n = 500
set.seed(20221224)
L = matrix(runif(p*q,max = 1,min = -1),nrow = p)# loading matrix
X = U\%*\%t(L) + matrix(rnorm(n*p,sd = 2),nrow = n)
Tuning_q(X)

%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.

%## The function is currently defined as
%function (x)
%{
%  }
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
}
