\name{qqpval}
\alias{qqpval}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
make a QQ plot of p-values.
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function makes QQ-plots for p-values, by default on a -log10 scale. It also shows pointwise-95\% confidence bounds for the order statistics of a Uniform(0,1) distribution.
}
\usage{
qqpval(p, cone = TRUE, log = TRUE, add = FALSE, col=1, pch=1,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{
  p-values
%%     ~~Describe \code{p} here~~
}
  \item{cone}{
  Logical. Whether or not to print 95\% confidence bounds for Uniform(0,1) order statistics
%%     ~~Describe \code{cone} here~~
}
  \item{log}{
  Logical. Whether or not to plot p-values on -log10 scale
%%     ~~Describe \code{log} here~~
}
  \item{add}{
  logical. whether or not to add to an existing plot
%%     ~~Describe \code{add} here~~
}
  \item{\dots}{
    other options to be passed to `plot'
%%     ~~Describe \code{\dots} here~~
}
  \item{col,pch}{
    color and point type. See \link{plot}
}
}
\author{
%%  ~~who you are~~
Arend Voorman
}
\examples{

p <- runif(1000)
qqpval(p)

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
}