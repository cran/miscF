mvt.mcmc <- function(X, prior.lower.v, prior.upper.v,
                     prior.Mu0=rep(0, ncol(X)),
                     prior.Sigma0=diag(10000, ncol(X)),
                     prior.p=ncol(X), prior.V=diag(1, ncol(X)),
                     initial.v=NULL, initial.Sigma=NULL,
                     nmcmc=10000, nburn=nmcmc/10, nthin=1, seed=1){
        
    if(!is.matrix(X) || ncol(X) < 2 || nrow(X) < 3){
        stop("The input observations of the 'mvt.mcmc' has to be a matrix
             of more than two rows and one column.")
    }
    if(prior.lower.v<0 || prior.upper.v<0 || prior.lower.v>prior.upper.v){
        stop("The bounds of the degrees of freedom of mvt have to be
             positive and the lower bound has to be smaller than the upper.")
    }
    if(!is.vector(prior.Mu0) || length(prior.Mu0)!=ncol(X)){
        stop("The mean vector of multivariate normal prior of the
             location of mvt has to be a vector of the same
             length as the number of columns of observations.")
    }
    if(!is.matrix(prior.Sigma0) || any(eigen(prior.Sigma0)$values<=0) || any(dim(prior.Sigma0)!=ncol(X))){
        stop("The variance matrix of multivariate normal
             prior of the location of mvt has to be a positive definite
             matrix with the number of rows and columns equal to
             the number of columns of observations.")
    }
    if(prior.p < ncol(X)){
        stop("The degrees of freedom of wishart prior of inverse of the
             scale matrix of mvt has be equal to or larger than the
             number of columns of observations.")
    }
    if(!is.matrix(prior.V) || any(eigen(prior.V)$values<=0) ||  any(dim(prior.V)!=ncol(X))){
        stop("The scale matrix of wishart prior of inverse of the scale
             matrix of mvt has to be a positive definite
             matrix with the number of rows and columns equal to
             the number of columns of observations.")
    }
    if(!is.null(initial.v) && (initial.v < prior.lower.v || initial.v > prior.upper.v)){
        stop("The initial value of the df of the mvt has to be within the bounds
             (prior.lower.v, prior.upper.v).")
    }
    if(!is.null(initial.Sigma) &&
           (!is.matrix(initial.Sigma) || any(eigen(initial.Sigma)$values<=0) || any(dim(initial.Sigma)!=ncol(X)))){
        stop("The initial value of the scale matrix of mvt has to be
             a positive definite matrix with the number of rows and columns
             equal to the number of columns of observations.")
    }
    
    
    n <- nrow(X)            
    d <- ncol(X)             
    
     
    #assign initial values obtained from ecme
    if(is.null(initial.v) || is.null(initial.Sigma)){
        ecme <- mvt.ecme(X, prior.lower.v, prior.upper.v)
        if(is.null(initial.v)){
            initial.v <- ecme$v
        }
        if(is.null(initial.Sigma)){
            initial.Sigma <- ecme$Sigma
        }
    }
    v <- initial.v
    Sigma <- initial.Sigma
    lambda <- rgamma(n, v/2, v/2)
    Sigma.inv <- solve(Sigma)
    
    datalist <- list(Y=X, n=n, k=d, 
                     Mu0=prior.Mu0, Tau0=solve(prior.Sigma0), H0=prior.V, p0=prior.p,
                     ndf1=prior.lower.v, ndf2=prior.upper.v)
    initlist <- list(list(nu=initial.v, Tau=solve(initial.Sigma)))
    parametersToSave <- c("Mu", "Sigma", "nu")
    
    cpath <- getwd()
    setwd(cpath)
    cat(
      "model {
            for( i in 1 : n ) {
                Y[i,1:k] ~ dmt(Mu[], Tau[,], nu)
            }
            Mu[1:k] ~ dmnorm(Mu0[], Tau0[,])
            Tau[1:k,1:k] ~ dwish(H0[,], p0)
            Sigma[1:k,1:k] <- inverse(Tau[,])
            nu ~ dunif(ndf1, ndf2)

       }",
        file="mvt_BUGS.txt")

    
    
    bugsfit <- BRugs::BRugsFit(modelFile="mvt_BUGS.txt", data=datalist, inits=initlist, numChains = 1, 
                    parametersToSave=parametersToSave,
                    nBurnin = nburn, nIter = nmcmc, nThin = nthin, coda = TRUE,
                    DIC = FALSE, working.directory = NULL, digits = 5, seed=seed)
    unlink("mvt_BUGS.txt")

    Mu <-  do.call(cbind, lapply(1:d, function(i) bugsfit[,paste0("Mu[", i, "]")][[1]]))
    Sigma <- array(0, dim=c(d, d, dim(bugsfit[[1]])[1]))
    for(i in 1:d){
      for(j in 1:d){        
        Sigma[i,j,] <- bugsfit[,paste0("Sigma[", i, ",", j, "]")][[1]]
      }
    }	  
    v <- bugsfit[,'nu'][[1]]
    
    list(Mu.save=Mu, Sigma.save=Sigma, v.save=v)

}
