\name{mvsn.mcmc}
\alias{mvsn.mcmc}
\title{Estimate Parameters of a Multivariate Skew Normal Distribution Using the MCMC 
}
\description{
  Use the MCMC to obtain estimate of parameters of a multivariate skew normal distribution.
}
\usage{
  mvsn.mcmc(Y, prior.Mu0=NULL, prior.Sigma0=NULL, 
            prior.muDelta0=NULL, prior.sigmaDelta0=NULL, 
            prior.H0=NULL, prior.P0=NULL,
            nmcmc=10000, nburn=nmcmc/10, nthin=1, seed=100)
}
\arguments{
  \item{Y}{a matrix of observations with one subject per row.}
  \item{prior.Mu0}{mean vector of multivariate normal prior of the
	parameter \eqn{\mu}. The default value is \code{NULL}. For the default,
        the value will be generated automatically.}
  \item{prior.Sigma0}{variance matrix of multivariate normal prior of
	the parameter \eqn{\mu}. The default value is \code{NULL}. For the default,
        the value will be generated automatically.}
  \item{prior.muDelta0}{mean vector of normal prior of the diagonal elements of parameter \eqn{D}. 
        The default value is \code{NULL}. For the default,
        the value will be generated automatically.}
  \item{prior.sigmaDelta0}{standard deviation vector of normal prior of the diagonal 
        elements of parameter \eqn{D}. 
        The default value is \code{NULL}. For the default,
        the value will be generated automatically.}
  \item{prior.H0}{the inverse of scale matrix of Wishart prior of the inverse of 
        parameter \eqn{\Sigma}. The default value is \code{NULL}. For the default,
        the value will be generated automatically.}
  \item{prior.P0}{the degrees of freedom of Wishart prior of the inverse of 
        parameter \eqn{\Sigma}. The default value is \code{NULL}. For the default,
        the value will be generated automatically.}
  \item{nmcmc}{number of iterations. The default value is 10000.}
  \item{nburn}{number of burn-in. The default value is \code{nmcmc}/10.}
  \item{nthin}{output every \code{nthin}-th sample. The default value is 1 (no thinning).}
  \item{seed}{random seed. The default value is 100.}
}
\value{
  \item{Mu}{a matrix of parameter \eqn{\mu} of the distribution, one row per iteration.}
  \item{Sigma}{a three dimensional array of parameter \eqn{\Sigma} of the
	distribution. Sigma[i,,] is the result from the i-th iteration.}
  \item{Delta}{a matrix of diagonal elements of parameter \eqn{D} of the distribution, one row per iteration.}
  \item{DIC}{DIC value.}
}
\details{
  This function estimates the parameters of a multivariate skew normal distribution as in Sahu et al. 2003 using the MCMC.
}
\references{
  Sahu, Sujit K., Dipak K. Dey, and Marcia D. Branco. (2003) 
  A new class of multivariate skew distributions with applications to Bayesian regression models.
  \emph{Canadian Journal of Statistics} \bold{vol. 31, no. 2} 129-150.
} 
\examples{
\dontrun{
  Mu <- rep(400, 2)
  Sigma <- diag(c(40, 40))
  D <- diag(c(-30, -30))
  Y <- rmvsn(n=1000, D, Mu, Sigma)
  mcmc <- mvsn.mcmc(Y)
}
}
\keyword{multivariate}



