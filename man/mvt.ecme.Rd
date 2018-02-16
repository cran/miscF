\name{mvt.ecme}
\alias{mvt.ecme}
\title{Estimate Parameters of a Multivariate t Distribution Using the
  ECME Algorithm}
\description{
  Use the Expectation/Conditional Maximization Either (ECME) algorithm
  to obtain estimate of parameters of a multivariate t distribution.
}
\usage{
  mvt.ecme(X, lower.v, upper.v, err=1e-4)
}
\arguments{
  \item{X}{a matrix of observations with one subject per row.}
  \item{lower.v}{lower bound of degrees of freedom (df).}
  \item{upper.v}{upper bound of df.}
  \item{err}{the iteration stops when consecutive difference in
             percentage of df reaches this bound. The default value is 1e-4.}
}
\value{
  \item{Mu}{estimate of location.}
  \item{Sigma}{estimate of scale matrix.}
  \item{v}{estimate of df.}
}
\details{
   They are number of forms of the generalization of the univariate student-t
   distribution to multivariate cases. This function adopts the widely
   used representation as a scale mixture of normal distributions.

   To obtain the estimate, the algorithm adopted is the
   Expectation/Conditional Maximization Either (ECME), which extends the
   Expectation/Conditional Maximization (ECM) algorithm by allowing
   CM-steps to maximize either the constrained expected complete-data
   log-likelihood, as with ECM, or the correspondingly constrained
   actual log-likelihood function.
 }
\references{
  Chuanhai Liu (1994)
  Statistical Analysis Using the Multivariate t Distribution
  \emph{Ph. D. Dissertation, Harvard University}
}
\examples{
  mu1 <- mu2 <- sigma12 <- sigma22 <- 100
  rho12 <- 0.7
  Sigma <- matrix(c(sigma12, rho12*sqrt(sigma12*sigma22),
                    rho12*sqrt(sigma12*sigma22), sigma22),
                  nrow=2)
  k <- 5
  N <- 100
  require(mvtnorm)
  X <- rmvt(N, sigma=Sigma, df=k, delta=c(mu1, mu2))

  result <- mvt.ecme(X, 3, 300)
  result$Mu
  result$Sigma
  result$v
}
\keyword{multivariate}



