mvst.mcmc <- function(Y, prior.Mu0=NULL, prior.Sigma0=NULL, 
                      prior.muDelta0=NULL, prior.sigmaDelta0=NULL, 
                      prior.H0=NULL, prior.P0=NULL,
        	      nmcmc=10000, nburn=nmcmc/10, nthin=1, seed=1){

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
 
    L.nu <- 4; U.nu <- 15 
    est.mvt <- mvt.ecme(Y, L.nu, U.nu, err=1e-4)
    nu.ini <- est.mvt$v
    init <- mvsn.mcmc(Y, nmcmc=1000)
    Delta.ini <- colMeans(init$Delta)
    Sigma.ini <- apply(init$Sigma, 2:3, mean)*(nu.ini-2)/nu.ini
    TauY.ini <- solve(Sigma.ini)

    datalist <- list(Y=Y, N=N, P=P, Mu0=Mu0, Tau0=Tau0, muDelta0=muDelta0, tauDelta0=tauDelta0, H0=H0, P0=P0)
    initlist <- list(list(Mu=Mu0, TauY=TauY.ini, Delta=Delta.ini, nu=nu.ini))
    parametersToSave <- c("Mu", "Delta", "Sigma", "nu")

    cpath <- getwd()
    setwd(cpath)
    cat(
    "model {
         for(i in 1:N){
             Y[i, 1:P] ~ dmnorm(MuY[i, 1:P], TauYW[i, 1:P, 1:P])
         }

         for(i in 1:N){
             for(j in 1:P){
	         MuY[i, j] <- Mu[j] + Dz[i, j]
             }
         }

         for(i in 1:N){
             for(j in 1:P){
                 for(k in 1:P){
	             TauYW[i, j, k] <- W[i]*TauY[j, k]
                 }
             }
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
                 Z[i, j] ~ dnorm(0, 1)I(0, )
             }
         }
	
         shape <- nu/2
         rate <- nu/2
         for(i in 1:N){
   	     W[i] ~ dgamma(shape, rate)
         }
	
         Mu[1:P] ~ dmnorm(Mu0[], Tau0[,]) 
    
         TauY[1:P,1:P] ~ dwish(H0[,], P0)
         Sigma[1:P,1:P] <- inverse(TauY[,])  

         nu ~ dexp(0.1)I(2, )
    }",
        file="mvst_BUGS.txt")

    BRugs::modelCheck("mvst_BUGS.txt")
    BRugs::modelData(BRugs::bugsData( datalist ))
    BRugs::modelCompile(numChains=1)
    BRugs::modelSetRN(seed)
    BRugs::modelInits(BRugs::bugsInits( initlist, numChains=1 ))
    BRugs::modelGenInits()    # to generate a starting value for the missing x
    BRugs::samplesSetThin(nthin)
    BRugs::modelUpdate(nburn)
    BRugs::dicSet()
    on.exit(BRugs::dicClear(), add = TRUE)
    BRugs::samplesSet(parametersToSave)
    BRugs::modelUpdate(nmcmc-nburn)

    bugsfit <- BRugs::buildMCMC("*")
    unlink("mvst_BUGS.txt")
	     
    Mu <-  do.call(cbind, lapply(1:P, function(i) bugsfit[,paste0("Mu[", i, "]")][[1]]))
    Delta <- do.call(cbind, lapply(1:P, function(i) bugsfit[,paste0("Delta[", i, "]")][[1]]))
    Sigma <- array(0, dim=c(dim(bugsfit[[1]])[1], P, P))
    for(i in 1:P){
        for(j in 1:P){        
	    Sigma[,i,j] <- bugsfit[,paste0("Sigma[", i, ",", j, "]")][[1]]
        }
    }	  
    nu <- bugsfit[,'nu'][[1]]
    DIC <- BRugs::dicStats()

    list(Mu=Mu, Delta=Delta, Sigma=Sigma, nu=nu, DIC=DIC[4,3])
    
}
