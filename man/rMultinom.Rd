\name{rMultinom}
\alias{rMultinom}
\title{Generate Random Samples from Different Multinomial
  Distributions with the Same Number of Classes
}
\description{
  Generate random samples from multinomial distributions with the same
  number of classes but different event probabilities.
}
\usage{
   rMultinom(p)
}
\arguments{
  \item{p}{matrix with each row specifying the probabilities for
	different classes of each sample.}
}
\value{
  Random samples from multinomial distributions.
}
\details{
  This function vectorizes the generation of random samples from
  different multinomial distributions by the inversion of CDF method.
}
\seealso{
  \code{\link[stats]{rmultinom}}
}
\examples{
  #Example 1: Generate 100 random samples from multinomial distributions
  #           with 3 classes and different event probabilities.
  p1 <- runif(100)
  p2  <- runif(100, 0, 1-p1)
  p3 <- 1-p1-p2 
  x <- rMultinom(p=cbind(p1, p2, p3))
  
}
\keyword{distribution}



