\name{mvn.bayes}
\alias{mvn.bayes}
\title{Estimate the Parameters of a Multivariate Normal Model by the
  Bayesian Methods}
\description{
  Estimate the parameters of a multivariate normal model under different
  priors.
}
\usage{
  mvn.bayes(X, nsim, prior=c("Jeffreys", "Conjugate"))
}
\arguments{
  \item{X}{a matrix of observations with one subject per row.}
  \item{nsim}{number of simulations.}
  \item{prior}{a character string specifying the prior distribution.
    It must be either "Jeffreys" or "Conjugate" and may be abbreviated.
	The default is "Jeffreys".}
}
\value{
  \item{Mu.save}{a matrix of mean vector of the model, one row per iteration.}
  \item{Sigma.save}{a three dimensional array of variance of the model.
	Sigma[,,i] is the result from the ith iteration.}
}
\details{
  Both the Jeffreys prior and normal-inverse-Wishart conjugate prior
  are available. The conjugate prior of variance covariance matrix is
  inverse-Wishart. To use a noninformative proper prior, the degree of
  freedom of Wishart prior was set as the number of dimensions and the
  scale matrix was chosen based on the unbiased estimate.
  The number of prior measurements was
  taken as one and the prior mean was set as its unbiased estimate.}
\note{When the number of dimensions is two, under Jeffreys prior,
  it is straightforward to obtain independent
  samples from the posteriors.
}
\references{
  Berger, J. O., Sun, D. (2008)
  Objective Priors for the Bivariate Normal Model. 
  \emph{The Annals of Statistics} \bold{36} 963-982.
 
  Gelman, A., Carlin, J. B., Stern, H. S., Rubin, D. B. (2003)
  Bayesian Data Analysis. 2nd ed. London: Chapman and Hall

  Sun, D., Berger, J. O. (2009)
  Objective Priors for the Multivariate Normal Model.
  In Bayesian Statistics 8, Ed. J. Bernardo, M. Bayarri, J. Berger,
  A. Dawid, D. Heckerman, A. Smith and M. West.
  Oxford: Oxford University Press.

} 
\examples{
  Sigma <- matrix(c(100, 0.99*sqrt(100*100),
                      0.99*sqrt(100*100), 100),
                      nrow=2)
  X <- mvrnorm(1000, c(100, 100), Sigma)
  result <- mvn.bayes(X, 10000)
  Mu <- colMeans(result$Mu.save)
  Sigma <- apply(result$Sigma.save, c(1,2), mean)
}
\keyword{multivariate}



