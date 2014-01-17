\name{diabetes}
\alias{diabetes}
\docType{data}
\title{
Blood and other measurements in diabetics.
}
\description{
The diabetes data frame has 442 rows and 3 columns. These are the data used in Efron et al "Least Angle Regression".
}
\usage{data(diabetes)}
\format{
  A data frame with 442 observations on the following 3 variables.
  \describe{
    \item{\code{x}}{a matrix with 10 columns}
    \item{\code{y}}{a numeric vector}
    \item{\code{x2}}{a matrix with 64 columns}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
The x matrix has been standardized to have variance 1 in each column and zero mean. The matrix x2 consists of x plus certain interactions.
}
\source{
Data can be found in the `lars' package.
}
\references{
Voorman, A, Shojaie, A, and Witten D. nference in high-dimensions with the penalized score test. In preparation.
}
\seealso{
\code{\link{lassoscore}}
}
\examples{
data(diabetes)
mod <- with(diabetes,lassoscore(y,x,lambda=0.02))
}
\keyword{datasets}
