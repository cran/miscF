\name{curve.polynomial.rjmcmc}
\alias{curve.polynomial.rjmcmc}
\title{Curve Fitting Using Piecewise Polynomials with Unknown Number and
  Location of Knots
}
\description{
  Fit a variety of curves by a sequence of piecewise polynomials.
  The number and location of knots are determined by the Reversible
  Jump MCMC method. 
}
\usage{
  curve.polynomial.rjmcmc(y, x, lambda, l, l0, c=0.4,
                          gamma.shape=1e-3, gamma.rate=1e-3,
                          maxit=10000, err=1e-8, verbose=TRUE)


}
\arguments{
  \item{y}{a vector of values of a response variable.}
  \item{x}{a vector of values of the corresponding explanatory variable.}
  \item{lambda}{the parameter of the Poisson prior for the number of knots.}
  \item{l}{the order of polynomials.}
  \item{l0}{the degree of continuity at the knots.}
  \item{c}{the parameter controlling the rate of changing dimension. The
           default value is 0.4.}
  \item{gamma.shape}{the parameter shape of the gamma prior for the error
                     precision. The default value is 1e-3.}
  \item{ gamma.rate}{the parameter rate of the gamma prior for the error
                     precision. The default value is 1e-3.}
  \item{maxit}{the maximum number of iterations. The default value is
               10000.}
  \item{err}{the iteration stops when consecutive difference in
             percentage of mean-squared error reaches this bound.
             The default value is 1e-8.}
  \item{verbose}{logical. If \code{TRUE}, then indicate
	the level of output after every 1000 iterations.
    The default is \code{TRUE}.}
}
\value{
   \item{knots.save}{a list of location of knots. One component per iteration.}
   \item{betas.save}{a list of estimates of regression parameters
                     \eqn{\beta}s. One component per iteration.}
   \item{fitted.save}{a matrix of fitted values. One column per iteration.}
   \item{sigma2.save}{a vector of estimate of error variance.
                      One component per iteration.}
}
\details{
   The method is based on Denison et al. (1998).  It can be used to fit
   both smooth and unsmooth curves.  When both \code{l0} and \code{l}
   are set to 3, it fits curves using cubic spline.
}
\note{
   The factor \eqn{\frac{1}{\sqrt{n}}} was added in the likelihood ratio
   to penalize the ratio for dimensionality as suggested in Dimatteo et
   al. (2001).
}
\references{
 Denison, D. G. T., Mallick, B. K., Smith, A. F. M.(1998) 
 Automatic Bayesian Curve Fitting
 \emph{JRSSB} \bold{vol. 60, no. 2} 333-350

 Dimatteo, I., Genovese, C. R., Kass, R. E.(2001) 
 Bayesian Curve-fitting with Free-knot Splines
 \emph{Biometrika} \bold{vol. 88, no. 4} 1055-1071
 
} 
\seealso{
  \code{\link[pspline]{sm.spline}}
}
\examples{
\dontrun{
   #Example 1: smooth curve
   #example 3.1. (b) in Denison et al.(1998)

   x <- runif(200)
   xx <- -2 + 4*x
   y.truth <- sin(2*xx) + 2*exp(-16*xx^2)
   y <- y.truth + rnorm(200, mean=0, sd=0.3)

   results <- curve.polynomial.rjmcmc(y, x, lambda=1, l=2, l0=1)

   plot(sort(x), y.truth[order(x)], type="l")
   lines(sort(x), rowMeans(results$fitted.save), col="red")

   #Example 2: unsmooth curve
   #blocks in Denison et al.(1998)
   tj <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)
   hj <- c(4, -5, 3, -4, 5, -4.2, 2.1, 4.3, -3.1, 2.1, -4.2)

   t <- seq(0, 1, len=2048)
   Ktmtj <- outer(t, tj, function(a,b) ifelse(a-b > 0, 1, -1))
   ft <- rowSums(Ktmtj \%*\% diag(hj))
   x <- t
   y <- ft + rnorm(2048, 0, 1)

   results <- curve.polynomial.rjmcmc(y, x, lambda=5, l=2, l0=1)

   plot(x, ft, type="s")
   lines(x, rowMeans(results$fitted.save), col="red")

   #Example 3: unsmooth curve
   #bumps in Denison et al.(1998)
   tj <- c(0.1, 0.13, 0.15, 0.23, 0.25, 0.4, 0.44, 0.65, 0.76, 0.78, 0.81)
   hj <- c(4, 5, 3, 4, 5, 4.2, 2.1, 4.3, 3.1, 5.1, 4.2)*10
   wj <- c(0.005, 0.005, 0.006, 0.01, 0.01, 0.03, 0.01, 0.01, 0.005, 0.008, 0.005)

   t <- seq(0, 1, len=2048)
   ft <- rowSums((1 + abs(outer(t, tj ,"-") \%*\% diag(1/wj)))^(-4) \%*\% diag(hj))
   y <- ft + rnorm(2048, 0, 1)

   results <- curve.polynomial.rjmcmc(y, t, lambda=5, l=2, l0=1)

   plot(t, ft, type="s")
   lines(t, rowMeans(results$fitted.save), col="red")
  }
}
\keyword{smooth}



