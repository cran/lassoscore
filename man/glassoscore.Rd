\name{glassoscore}
\alias{glassoscore}
\alias{mbscore}

\title{
Penalized score test, for the graphical lasso.
}
\description{
Test whether an element of the precision matrix is zero, using the graphical lasso to approximate the parameters in remainder of the precision matrix.
}
\usage{
glassoscore(x, lambda, subset=NULL, penalize.diagonal=FALSE, tol=1e-8)
mbscore(x, lambda, subset=NULL, tol=1e-8,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  data matrix. Unlike \link{glasso}, this function requires the original data, not just the covariance matrix.
}
  \item{lambda}{
  a non-negative tuning parameter
}
\item{subset}{
  An ncol(x) by ncol(x) logical matrix, giving a subset of edges to test.
}
\item{penalize.diagonal}{
logical. Whether or not to penalize the diagonal in the graphical lasso. Defaults 
to FALSE.
}
  \item{tol}{
  convergence tolerance for \link{glasso} or \link{glmnet}
  }
  \item{...}{for mbscore, additional arguments to be passed to \link{lassoscore}
}
}
\details{
  This function tests for pairwise association between features, using the graphical lasso (glassoscore) or neighborhood selection (mbscore). Tests are based on the penalized score statistic \eqn{T[\lambda]}{T_\lambda}, described in Voorman \emph{et al} (2014). Note that a feature is non-zero in the (graphical) lasso solution if and only if \deqn{| T[\lambda] | > \lambda \sqrt{n},}{| T_\lambda | > \lambda \sqrt n,}
  where \eqn{T[\lambda]}{T_\lambda} is penalized the score statistic.
  
Calculating the variance of \eqn{T[\lambda]}{T_\lambda} can be computationally expensive for glassoscore. If there are \eqn{q} non-zero parameters in the graphical lasso solution, it will (roughly) require construction, and inversion, of a \eqn{q \times q}{q x q} matrix for each of the \eqn{q} non-zero parameters. That is, complexity is roughly \eqn{q^4}.

For mbscore, the results are typically not symmetric. For instance, p.sand[-i,i] contains the p-values produced by \code{lassoscore(x[,i],x[,-i],lambda)}, i.e. using x[,i] as the outcome variable, and thus p.sand[i,-i] contains p-values associated with feature i when used as the a predictor variable.
}
\value{
for an object of class either `glassoscore' or `mbscore', containing
  \item{scores}{the penalized score statistics}
  \item{scorevar.model}{the variance of the score statistics, estimated using a model-based variance estimate}
  \item{scorevar.sand}{the variance of the score statistcs, using a conservative variance estimate}
  \item{p.model}{p-value, using the model-based variance}
  \item{p.sand}{p-value, using the sandwich variance}
  \item{beta}{for mbscore, the beta[-i,i] contains the coefficients from lasso regression of x[,i] on x[,-i].}
In addition, glassoscore contains the output from `glasso' applied to x.
}
\references{
Jerome Friedman, Trevor Hastie and Robert Tibshirani (2007). Sparse inverse covariance estimation with the lasso. Biostatistics 2007. \url{http://www-stat.stanford.edu/~tibs/ftp/graph.pdf}

N. Meinshausen and P. Buhlmann. High-dimensional graphs and variable selection with the lasso. Annals of Statistics, 34(3):1436-1462, 2006.
}
\author{
Arend Voorman
}
\seealso{
\link{lassoscore}, \link{glasso}
}
\examples{
set.seed(100)

x<-matrix(rnorm(50*20),ncol=20)

gl <- glassoscore(x,0.2)
mb <- mbscore(x,0.2)

par(mfrow=c(1,2))
plot(gl)
plot(mb)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.

