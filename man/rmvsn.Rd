\name{rmvsn}
\alias{rmvsn}
\title{Generate Random Samples from a Multivariate Skew Normal Distribution
}
\description{
  Generate random samples from a multivariate skew normal distribution.
}
\usage{
   rmvsn(n, D, Mu, Sigma)
}
\arguments{
  \item{n}{number of samples.}
  \item{D}{parameter \eqn{D} of the distribution.}
  \item{Mu}{parameter \eqn{\mu} of the distribution.}
  \item{Sigma}{parameter \eqn{\Sigma} of the distribution.}
}
\value{
  Random samples from the multivariate skew normal distribution.
}
\details{
  This function generates random samples using the methods in Sahu et al. 2003.
}
\references{
  Sahu, Sujit K., Dipak K. Dey, and Marcia D. Branco. (2003) 
  A new class of multivariate skew distributions with applications to Bayesian regression models.
  \emph{Canadian Journal of Statistics} \bold{vol. 31, no. 2} 129-150.
} 
\examples{
  Mu <- rep(400, 2)
  Sigma <- diag(c(40, 40))
  D <- diag(c(-30, -30))
  Y <- rmvsn(n=1000, D, Mu, Sigma)
}
\keyword{distribution}



