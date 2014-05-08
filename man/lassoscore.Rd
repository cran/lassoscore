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
Test for the association between y and each column of X, adjusted for the other columns using a lasso regression, as described in Voorman et al (2014).
}
\usage{
lassoscore(y,X, lambda=0, family=c("gaussian","binomial","poisson"), 
    tol = .Machine$double.eps, maxit=1000, 
    resvar = NULL, verbose=FALSE, subset = NULL)
}
\arguments{
  \item{y}{
	outcome variable
}
  \item{X}{
	matrix of predictors
}
  \item{lambda}{
	tuning parameter value (see details)
}
\item{family}{
  The family, for the likelihood.
}
  \item{tol,maxit}{
	convergence tolerance and maximum number of iterations in \code{\link{glmnet}}
}
  \item{resvar}{
  value for the residual variance, for "gaussian" family. If not specified, the residual variance from lasso regression on all features is used (see details).
}
  \item{verbose}{
	 whether or not to print progress bars (defaults to FALSE)
}
  \item{subset}{
	 a subset of columns to test
}
}
\details{
For each column of X, denoted by x*, this function computes the score statistic 

\deqn{T[\lambda] = x*^T(y-\hat y)/\sqrt{n} }{T_\lambda = x*^T(y- yhat)/\sqrt n,} 

where \eqn{yhat} are the fitted values from lasso regression of y on X[,-x*] (see Note 2).

The variance of the score statistic is estimated in 4 ways: 

(i) a model-based estimate 

(ii) a sandwich varaince

(iii/iv) conservative versions of (i) and (ii), which do not depend on the selected model

Note 1: in lasso regression of y on X, the coefficient of x* is non-zero \bold{if and only if} 

\deqn{| T[\lambda] | > \lambda \sqrt{n}}{| T_\lambda | > \lambda \sqrt n}

Note 2: For lasso regression of y on X, we minimize -l(b) + lambda*||b||_1 over vectors b, where l(b) is either RSS/(2n) (for the "gaussian" family), or the log-likelihood for a generalized linear model. See the details of \link{glmnet} for more information.

Note 3:Each feature x is rescaled to have mean zero and x^Tx/n = 1, y is centered, but not rescaled.

}
\value{
Object of class `lassoscore', which is an R `list', with elements:
  \item{fit}{Elements of the fitted lasso regression of y on X (see \code{\link{glmnet}} for details.)} 
  \item{scores}{the score statistics}
  \item{resvar}{the value used for the residual variance}
  
  \item{scorevar.model}{the variance of the score statistics, estimated using a model-based approximation}
  \item{scorevar.sand}{the variance of the score statistcs, using an model-agnostic, or sandwich formula}
  \item{scorevar.model.cons,scorevar.sand.cons}{conservative versions of the variances}
  
 \item{p.model}{p-value, using a model-based variance}
 \item{p.sand}{p-value, using sandwich variance}
 \item{p.model.cons,p.sand.cons}{p-value, using conservative variance formulas}
}
\references{
Voorman, A, Shojaie, A, and Witten D. Inference in high dimensions with the penalized score test. \url{http://arxiv.org/abs/1401.2678}.
}
\author{
Arend Voorman \email{voorma@uw.edu}
}

\seealso{
\code{\link{glassoscore}}, \code{\link{qqpval}}
}
\examples{
#Simulation from Voorman et al (2014)
set.seed(20)
n <- 300
p <- 100
q <- 10

set.seed(20)
beta <- numeric(p)
beta[sample(p,q)] <- 0.4

Sigma <- forceSymmetric(t(0.5^outer(1:p,1:p,"-")))
cSigma <- chol(Sigma)

x <- scale(replicate(p,rnorm(n))\%*\%cSigma)
y <- rnorm(n,x\%*\%beta,1)

mod <- lassoscore(y,x,0.02)
summary(mod)
plot(mod,type="all")

#test only features 10:20:
mod0 <- lassoscore(y,x,0.02, subset = 10:20)

######## Diabetes data set:
#Test features in the diabetes data set, using 2 different values of `lambda', 
#and compare results:
resvar <- with(lm(y~x,data=diabetes), sum(residuals^2)/df.residual)

mod2 <- with(diabetes,lassoscore(y,x,lambda=4,resvar=resvar))
mod3 <- with(diabetes,lassoscore(y,x,lambda=0.5,resvar=resvar))
data.frame(
  "variable"=colnames(diabetes$x),
  "lambda_4"=format(mod2$p.model,digits=2),
  "lambda_0.5"=format(mod3$p.model,digits=2))
}
