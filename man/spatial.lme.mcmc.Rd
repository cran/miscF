\name{spatail.lme.mcmc}
\alias{spatial.lme.mcmc}
\title{Spatial Modeling by a Bayesian Hierarchical
       Linear Mixed-effects Model 
}
\description{
  A linear mixed-effects model that combines unstructured
  variance/covariance matrix for
  inter-regional (long-range) correlations and 
  an exchangeable correlation structure for intra-regional (short-range)
  correlations. Estimation is performed using the Gibbs sampling. 
}
\usage{
  spatial.lme.mcmc(spatialMat, nlr, nsweep, verbose=TRUE)
}
\arguments{
    \item{spatialMat}{a matrix of observations with one subject per column
      and one location per row. For each subject, the observations
	  are arranged one location after another.}
    \item{nlr}{a vector of number of locations within each region. One
	  component per region.}
	\item{nsweep}{the number of iterations.}
	\item{verbose}{logical. If \code{TRUE}, then indicate
	the level of output after every 1000 iterations.
    The default is \code{TRUE}.}

}
\value{
   \item{mu.save}{a matrix of means at every location. One column per iteration.}
   \item{sigma2.save}{a matrix of intra-regional variances. One column
	 per iteration.}
   \item{lambda2.save}{a matrix of variances of locations' means within regions. One column per iteration.}
   \item{Gamma.save}{a matrix of inter-regional variance/covariance
	 matrix. One column per iteration and within each column the
	 elements of the variance/covariance matrix are arranged column-wise.}
}
\details{
   The function was proposed to study the fMRI data.
   The original MATLAB code written by DuBois Bowman and Brian Caffo can
   be found at: \url{http://www.biostat.jhsph.edu/~bcaffo/downloads/clusterBayes.m}. 
   Instead of stacking the data from
   two conditions, the R version fits the model for one condition and
   the user needs to use the function multiple times for separate conditions.

   The initial values are obtained based on sample moments. The
   hyper-parameters for the prior distributions of the intra-regional
   variances and variances of locations' means are set up in the way that
   the mean is equal to the sample mean and the variance is large.
}
\note{
  There seemed to be no easy way to use \code{\link[lme4]{lmer}} or
  \code{\link[nlme]{lme}} to fit
  the variance/covariance structure in this model and \emph{SAS proc
  mixed} failed for certain cases.
  
}
\references{
  Brian Caffo, DuBois Bowman, Lynn Eberly and Susan Spear Bassett (2009)
  A Markov Chain Monte Carlo Based Analysis of a
  Multilevel Model for Functional MRI Data
  \emph{Handbook of Markov Chain Monte Carlo}

  F. DuBois Bowman, Brian Caffo, Susan Spear Bassett, and Clinton Kilts (2008)
  A Bayesian Hierarchical Framework for Spatial Modeling of fMRI
  Data
  \emph{Neuroimage} \bold{vol. 39, no. 1} 146-156
} 
\examples{
\dontrun{
    #simulate the data
    ns=100; nr=2; nlr <- c(20, 20)
    mu0 <- c(0, 0)
    sigma2 <- c(1., 1.)
    Gamma <- matrix(c(3, 0, 0, 3), nrow=2)

    sample <- matrix(0, nrow=sum(nlr), ncol=ns)
    for(i in 1:ns){
        alpha <- mvrnorm(1, rep(0, nr), Gamma)
        sampleR <- NULL
        for(g in 1:nr){
            beta <- rnorm(nlr[g], mean=alpha[g] + mu0[g], sd=sqrt(sigma2[g]))
            sampleR <- c(sampleR, beta)
        }
        sample[,i] <- sampleR
    }

    #run mcmc
    mcmc.result <- spatial.lme.mcmc(sample, nlr, 10000)

    #check the results
    Gamma <- mcmc.result$Gamma.save
    sigma2 <- mcmc.result$sigma2.save
    mu <- mcmc.result$mu.save
    matrix(rowMeans(Gamma), nr, nr)
    apply(sigma2, 1, function(x) quantile(x, prob=c(0.025, 0.5, 0.975)))
    summary(rowMeans(mu[1:nlr[1],]))
    summary(rowMeans(mu[(nlr[1]+1):sum(nlr),]))
    
}
}
\keyword{spatial}



