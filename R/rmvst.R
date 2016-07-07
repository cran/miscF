rmvst <- function(n, D, Mu, Sigma, nu){
    if(!is.matrix(D) || nrow(D) != ncol(D))
        stop("'D' has to be a square matrix")
    if(any(D[upper.tri(D)] != 0) | any(D[lower.tri(D)] != 0))
        stop("'D' has to be a diagnonal matrix")
    if(!is.vector(Mu) || length(Mu)!=nrow(D))
        stop("'Mu' has to a vector with the length equal to the number of rows of 'D'.")
    if(! isSymmetric(Sigma))
        stop("'Sigma' has to be a symmetric matrix.")
    if(! all(eigen(Sigma)$values > 0))
        stop("'Sigma' has to be positive definite.")
    if(nrow(Sigma) != nrow(D))
        stop("The dimension of 'Sigma' does not match that of 'D'.")
    if(nu <= 2)
        stop("'nu' has to be larger than two.")

    p <- nrow(D)
    Y <- matrix(0, n, p)  

    for(i in 1:n){
        z <- rtnorm(p, 0, 1, left=0)
        w <- rgamma(1, nu/2, nu/2)
        Y[i,] <- mvrnorm(1, Mu+D%*%z, Sigma/w)
    }

    Y
}
