\name{uvnm.rjmcmc}
\alias{uvnm.rjmcmc}
\title{Univariate Normal Mixture (UVNM) Model with Unknown Number of Components 
}
\description{
  Estimate the parameters of an univariate normal mixture model including
  the number of components using the Reversible Jump MCMC method. It can
  be used for density estimation and/or classification.
}
\usage{
  uvnm.rjmcmc(y, nsweep, kmax, k, w, mu, sigma2, Z,
              delta=1, xi=NULL, kappa=NULL, alpha=2,
              beta=NULL, g=0.2, h=NULL, verbose=TRUE)
}
\arguments{
  \item{y}{a vector of observations.}
  \item{nsweep}{number of sweeps. One sweep has six moves to update
	parameters including the number of components.}
  \item{kmax}{the maximum number of components.}
  \item{k}{initial value of number of components.}
  \item{w}{initial values of weights of different components.}
  \item{mu}{initial values of means of different components.}
  \item{sigma2}{initial values of variances of different components.}
  \item{Z}{initial values of allocations of all observations.}
  \item{delta}{the parameter of the prior distribution of \code{w},
	the symmetric Dirichlet distribution. The default value is one.}
  \item{xi}{the mean of the prior distribution of means of
	components. By taking the default value \code{NULL}, it is
	set as the median of the data \code{y} internally.}
  \item{kappa}{the precision of the prior distribution of means of
	components. By taking the default value \code{NULL}, 
	it is set as 1/R^2 internally, where R is the range of \code{y}.}
  \item{alpha}{the parameter shape of the prior distribution of
	precision of components. The default values is 2.}
  \item{beta}{the parameter rate of the prior distribution of
	precision of components. By taking the default value 
	\code{NULL}, it is set as
    a number generated randomly from a gamma distribution with shape
	\code{g} and rate \code{h}.}
  \item{g}{the parameter shape of the gamma distribution for
	\code{beta}. The default value is 0.2.}
  \item{h}{the parameter rate of the gamma distribution for
	\code{beta}. By taking the default value \code{NULL}, it is
	set as 10/R^2 internally, where R is the range of \code{y}.}
  \item{verbose}{logical. If \code{TRUE}, then indicate
	the level of output after every 1000 sweeps.
    The default is \code{TRUE}.}

}
\value{
  \item{k.save}{a vector of number of components.}
  \item{w.save}{weights of the UVNM, one component of the list per sweep.}
  \item{mu.save}{means of the UVNM, one component of the list per sweep.}
  \item{sigma2.save}{variances of the UVNM, one component of the list per sweep.}
  \item{Z.save}{a matrix of allocation, one column per sweep.}
}
\details{
    Estimate the parameters of a univariate normal mixture model with flexible
    number of components using the Reversible Jump MCMC method in
	Richardson and Green (1997).
}
\note{
  The error in equation (12) of Richardson and Green (1997)
  was corrected based on Richardson and Green (1998).
 
}
\references{
 Sylvia Richardson and Peter J. Green (1997) 
 On Bayesian Analysis of Mixtures with an Unknown Number of Components
 \emph{JRSSB} \bold{vol. 59, no. 4} 731-792

 Sylvia Richardson and Peter J. Green (1998) 
 Corrigendum: On Bayesian Analysis of Mixtures with an Unknown Number of Components
 \emph{JRSSB} \bold{vol. 60, no. 3} 661
 
} 
\seealso{
  \code{\link[mixAK]{NMixMCMC}}
}
\examples{

\dontrun{
    require(mixAK)
    data(Acidity)
    y <- Acidity
    w <- c(0.50, 0.17, 0.33)
    mu <- c(4, 5, 6)
    sigma2 <- c(0.08, 0.10, 0.14)
    Z <- do.call(cbind, lapply(1:3, function(i)
                                    w[i]*dnorm(y, mu[i], sqrt(sigma2[i]))))
    Z <- apply(Z, 1, function(x) which(x==max(x))[1])
     
    result <- uvnm.rjmcmc(y, nsweep=200000, kmax=30, k=3,
                          w, mu, sigma2, Z)

    ksave <- result$k.save
    round(table(ksave[-(1:100000)])/100000,4)

    #conditional density estimation
    focus.k <- 3
    pick.k <- which(ksave==focus.k)
    w <- unlist(result$w.save[pick.k])
    mu <- unlist(result$mu.save[pick.k])
    sigma2 <- unlist(result$sigma2.save[pick.k])
    den.estimate <- rep(w, each=length(y)) *
                    dnorm(rep(y, length(w)), mean=rep(mu, each=length(y)),
                          sd=rep(sqrt(sigma2), each=length(y)))
    den.estimate <- rowMeans(matrix(den.estimate, nrow=length(y)))*focus.k

    #within-sample classification
    class <- apply(result$Z.save[,pick.k], 1,
                                function(x) c(sum(x==1), sum(x==2), sum(x==3)))
    class <- max.col(t(class))

    #visualize the results
    hist(y, freq=FALSE, breaks=20, axes=FALSE, ylim=c(-0.3, 1),
         main="Density Estimation and Classification", ylab="")
    axis(2, at=c(-(3:1)/10, seq(0,10,2)/10), labels=c(3:1, seq(0,10,2)/10),
         font=2)
    lines(sort(y), den.estimate[order(y)], col="red")
    for(i in 1:3){
        points(y[class==i], rep(-i/10, length(y[class==i])), col=i, pch=i)
    }
    mtext("Density", 2, at=0.5, line=2)
    mtext("Classification", 2, at=-0.15, line=2)
}
}
\keyword{distribution}



