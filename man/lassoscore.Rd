\name{lassoscore}
\alias{lassoscore}
\alias{lassoscore-package}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
Lasso penalized score test
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
Test for the association between y and each column of X, adjusted for the other columns using a lasso regression, as described in Voorman et al (2013).
}
\usage{
lassoscore(y, X, lambda = 0, tol = .Machine$double.eps, 
  maxit = 1000, verbose = FALSE, subset = NULL, resvar = NULL)
}
\arguments{
  \item{y}{
	outcome variable (continuous)
}
  \item{X}{
	matrix of predictors
}
  \item{lambda}{
	tuning parameter value (see details)
}
  \item{tol}{
	convergence tolerance (see details in \code{\link{lasso}})
}
  \item{maxit}{
	maximum number of iterations for coordinate descent lasso algorithm (see details in \code{\link{lasso}})
}
  \item{verbose}{
	 whether or not to print progress bars (defaults to FALSE)
}
  \item{subset}{
	 a subset of variables to test
}
  \item{resvar}{
	value for the residual variance. If not specified, the residual variance from lasso regression on all features is used (see details).
}
}
\details{
For each column of X, denoted by x*, this function computes the score statistic 

T = (x*^T(y-yhat))/sqrt(n), 

where `yhat' are the fitted values from lasso regression of y on X[,-x*] (see Note 2).

The variance of the score statistic is estimated in 3 ways: 

(i) a conservative approximation, 

(ii) asymtotic variance, assuming non-random X,  and 

(iii) a sandwich, (or model-agnostic) formulation.

Note 1: in lasso regression of y on X, the coefficient of x* is non-zero if and only if |T| > sqrt(n)*lambda

Note 2: For lasso regression of y on X, we minimize ||y - Xb||_2^2/(2n) + lambda*||b||_1 over vectors b.

Note 3:Each feature x is rescaled to have mean zero and x^Tx/n = 1, y is centered, but not rescaled.

}
\value{
Object of class `lassoscore', which is an R `list', with elements:
  \item{fit}{Elements of the fitted lasso regression of y on X (see \code{\link{lasso}} for details.)} 
  \item{scores}{the score statistics}
  \item{resvar}{the value used for the residual variance}
  \item{scorevar.cons}{the variance of the score statistics, estimated using a conservative approximation}
  \item{scorevar.asm}{the variance of the score statistcs, using an asymtotic estimate assuming non-random X}
  \item{scorevar.sand}{the variance of the score statistics, using model-agnostic (sandwich) estimate}
  \item{p.cons}{p-value, using a conservative approximation}
  \item{p.asm}{p-value, assuming fixed design, and constant variance.}
  \item{p.sand}{p-value, using model-agnostic variance (sandwich) calculations.}
}
\references{
Voorman, A, Shojaie, A, and Witten D. Inference in high dimensions with the penalized score test. submitted. \url{http://arxiv.org/abs/1401.2678}.
}
\author{
Arend Voorman
}

\seealso{
\code{\link{lasso}}, \code{\link{qqpval}}
}
\examples{
#Simulation from in Voorman et al (2013)
set.seed(20)
n <- 300
p <- 100
q <- 8

Sigma <- matrix(0,p-q,p-q)
diag(Sigma) <- 1
for(i in 1:(p-q-2)){
  Sigma[i,i+1] <- Sigma[i+1,i] <- 0.5
  Sigma[i,i+2] <- Sigma[i+2,i] <- 0.25
}
Sigma[p-q-1,p-q] <- Sigma[p-q,p-q-1] <- 0.5

X <- replicate(p,rnorm(n))
X[,-(1:q)] <- as.matrix(X[,-(1:q)]\%*\%chol(Sigma))
X <- scale(X)*sqrt(n/(n-1))

beta <- c(rep(.15,q), rep(0,p-q))
y <- scale(rnorm(n, X\%*\%beta),scale=FALSE)

mod <- lassoscore(y,X,0.02)
summary(mod)
plot(mod,type="all")

#test only features 10:20:
mod0 <- lassoscore(y,X,0.02, subset = 10:20)

######## Diabetes data set:
#Test features in the diabetes data set, using 2 different values of `lambda', 
#and compare results:
mod2 <- with(diabetes,lassoscore(y,x,lambda=7))
mod3 <- with(diabetes,lassoscore(y,x,lambda=1))

data.frame(
  "variable"=colnames(diabetes$x),
  "lambda_4"=format(mod2$p.asm,digits=3),
  "lambda_0.2"=format(mod3$p.asm,digits=3))

}
