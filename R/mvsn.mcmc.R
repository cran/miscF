mvsn.mcmc <- function(Y, prior.Mu0=NULL, prior.Sigma0=NULL, 
                      prior.muDelta0=NULL, prior.sigmaDelta0=NULL, 
                      prior.H0=NULL, prior.P0=NULL,
                      nmcmc=10000, nburn=nmcmc/10, nthin=1, seed=100){

    if(!is.matrix(Y) || ncol(Y) < 2 || nrow(Y) < 3){
        stop("The input observations of the 'mvsn.mcmc' has to be a matrix
             of more than two rows and one column.")
    }
   
    if(!is.null(prior.Mu0) && (!is.vector(prior.Mu0) || length(prior.Mu0)!=ncol(Y))){
        stop("The mean vector of multivariate normal prior of the
             parameter 'mu' of mvsn has to be a vector of the same
             length as the number of columns of observations.")
    }

    if(!is.null(prior.Sigma0) && (!is.matrix(prior.Sigma0) || any(eigen(prior.Sigma0)$values<=0) || any(dim(prior.Sigma0)!=ncol(Y)))){
        stop("The variance matrix of multivariate normal
             prior of the parameter 'mu' of mvsn has to be a positive definite
             matrix with the number of rows and columns equal to
             the number of columns of observations.")
    }
  
    if(!is.null(prior.muDelta0) && (!is.vector(prior.muDelta0) || length(prior.sigmaDelta0)!=ncol(Y))){
        stop("The vector of the means of normal prior
             of parameter 'D' of mvsn has to be a vector of the same
             length as the number of columns of observations.")
    }

    if(!is.null(prior.sigmaDelta0) && (!is.vector(prior.sigmaDelta0) || length(prior.sigmaDelta0)!=ncol(Y) || any(prior.sigmaDelta0 <= 0))){
        stop("The vector of the standard deviations of normal prior
             of parameter 'D' of mvsn has to be a vector of the same
             length as the number of columns of observations and 
             its components have to be positive.")
    }

    if(!is.null(prior.H0) && (!is.matrix(prior.H0) || any(eigen(prior.H0)$values<=0) || any(dim(prior.H0)!=ncol(Y)))){
        stop("The inverse of scale matrix of Wishart prior of the inverse of 
             parameter 'Sigma' of mvsn has to be a positive definite
             matrix with the number of rows and columns equal to
             the number of columns of observations.")
    }
 
    	
    if(!is.null(prior.P0) && (!is.vector(prior.P0) || length(prior.H0)!=1 || prior.H0 < ncol(Y))){
        stop("The degrees of freedom of Wishart prior of the inverse of 
             parameter 'Sigma' of mvsn has to be positive scalar and no less than
             the number of columns of observations.")
    }

     

    P <- ncol(Y)
    N <- nrow(Y)

    if(is.null(prior.Mu0)){
        Mu0 <- apply(Y, 2, median)
    }
    else{
	Mu0 <- prior.Mu0
    }

    if(is.null(prior.Sigma0)){
	prior.Sigma0 <- 100*var(Y)  
    }
    Tau0 <- solve(prior.Sigma0)
    
    if(is.null(prior.muDelta0)){
	prior.muDelta0 <- rep(0, P)
    }
    muDelta0 <- prior.muDelta0

    if(is.null(prior.sigmaDelta0)){
	prior.sigmaDelta0 <- sqrt(diag(var(Y)))
    }
    tauDelta0 <- 1/(prior.sigmaDelta0^2)
    
    if(is.null(prior.P0)){ 
        P0 <- P
    }
    else{
	P0 <- prior.P0
    }

    if(is.null(prior.H0)){
	H0 <- P*var(Y)
    }
    else{
	H0 <- prior.H0
    }
    
    Delta.ini <- mvrnorm(1, rep(0, P), diag(diag(var(Y)/(1-2/pi))))

    data=list("Y", "N", "P", "muDelta0", "tauDelta0", "Mu0", "Tau0", "H0", "P0")
    inits = function(){list(Delta=Delta.ini)}
    
    cpath <- getwd()
    setwd(cpath)
    cat(
    "model {
    
        for(i in 1:N){
    	    Y[i,1:P] ~ dmnorm (MuY[i, 1:P], TauY[1:P, 1:P])
            MuY[i, 1:P] <- Mu + Dz[i, 1:P]
	    }

        for(i in 1:N){
            for(j in 1:P){
                Dz[i,j] <- Delta[j]*Z[i,j]
            }
        }

        for(j in 1:P){
	    Delta[j] ~ dnorm(muDelta0[j], tauDelta0[j])
        }	 

        for(i in 1:N){
            for(j in 1:P){
                Z[i, j] ~ dnorm(0, 1)T(0, )
            }
        }
    
        Mu[1:P] ~ dmnorm(Mu0[1:P], Tau0[1:P, 1:P]) 
	
        TauY[1:P , 1:P] ~ dwish(H0[1:P, 1:P], P0)  
        Sigma[1:P, 1:P] <- inverse(TauY)
    }",
        file="mvsn.jags")

    set.seed(seed)	
    jagsfit <- jags(data, inits,
                    parameters.to.save=c("Mu", "Delta", "Sigma"), 
                    model.file="mvsn.jags",
                    n.chains = 1, n.iter=nmcmc, n.burnin=nburn, n.thin=nthin)
    unlink("mvsn.jags")

    Mu <- jagsfit$BUGSoutput$sims.list$Mu
    Sigma <- jagsfit$BUGSoutput$sims.list$Sigma
    Delta <- jagsfit$BUGSoutput$sims.list$Delta
    DIC=jagsfit$BUGSoutput$DIC
    
    list(Mu=Mu, Delta=Delta, Sigma=Sigma, DIC=DIC)

}
