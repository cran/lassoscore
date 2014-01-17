\name{lasso}
\alias{lasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
fit lasso regression
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This is a simple lasso implementation, meant to be used internally for lasso score test. This function minimizes ||y-Xb||_2^2/(2n) + lambda*||b||_1, over vectors b, using coordinate descent. For better, and more general, lasso implementations, see e.g. the `plus' or 'glmnet' packages.
}
\usage{
lasso(y, X, lambda, beta = NULL, maxit = 1000, tol = 2e-16, rescale=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
  an outcome variable
%%     ~~Describe \code{y} here~~
}
  \item{X}{
  a matrix of predictors
%%     ~~Describe \code{X} here~~
}
  \item{lambda}{
  The tuning parameter
%%     ~~Describe \code{lambda} here~~
}
  \item{beta}{
  initial coefficient values (optional)
%%     ~~Describe \code{beta} here~~
}
  \item{maxit}{
%%     ~~Describe \code{maxit} here~~
  maximum number of iterations the lasso algorithm (see details).
}
  \item{tol}{
%%     ~~Describe \code{tol} here~~
convergence tolerance (see details)
}
\item{rescale}{
whether or not to rescale x to have mean-zero and variance 1.
}
}
\details{
This function implements a simple coordinate descent algorithm for lasso regression. At each iteration, the algorithm cycles through the features j=1...d, minimizing ||y-Xb||_2^2/(2n) + lambda*||b||_1 with respect to b[j]. The algorithm terminates if the maximum change in the coefficient estimate is smaller than `tol' or the number of iterations reaches `maxit'.
}
\value{
%%  ~Describe the value returned
An R `list', with elements:
  \item{r}{residuals} 
	\item{X}{predictor matrix, as a numeric vector} 
	\item{y}{outcome} 
	\item{beta}{coefficient}
	\item{p}{dimension (ncol(X))} 
	\item{n}{sample size} 
	\item{lambda}{tuning parameter value}
	\item{maxit}{maximum number of iterations}
	\item{eps}{convergence tolerance}
}
\author{
Arend Voorman
}
\seealso{
\code{\link{lassoscore}}
}
