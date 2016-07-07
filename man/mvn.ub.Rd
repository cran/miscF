\name{mvn.ub}
\alias{mvn.ub}
\title{Unbiased Estimate of Parameters of a Multivariate Normal Distribution}
\description{
  Obtain the Unbiased Estimate of Parameters of a Multivariate Normal Distribution.
}
\usage{
  mvn.ub(X)
}
\arguments{
  \item{X}{a matrix of observations with one subject per row.}
}
\value{
  \item{hat.Mu}{unbiased estimate of mean.}
  \item{hat.Sigma}{unbiased estimate of variance.}
}
\examples{
  Sigma <- matrix(c(100, 0.99*sqrt(100*100),
                      0.99*sqrt(100*100), 100),
                      nrow=2)
  X <- mvrnorm(1000, c(100, 100), Sigma)
  result <- mvn.ub(X)
}
\keyword{multivariate}



