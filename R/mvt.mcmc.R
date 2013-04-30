#log likelihood of df v
logv <- function(v, n, lambda){
	vh <- v/2
    n*vh*log(vh) - n*lgamma(vh) + (vh-1)*sum(log(lambda)) - vh*sum(lambda)
}

#slice sampling of v
v.slice <- function(v0, n, lambda, logv, w, m, lower, upper, logv0=NULL)
{
	uni.slice.evals <- 0
	
	if (is.null(logv0)){
        uni.slice.evals <<- uni.slice.evals + 1
		logv0 <- logv(v0, n, lambda)
	}
	
	logy <- logv0 - rexp(1)
	
	u <- runif(1,0,w)
	L <- v0 - u
	R <- v0 + (w-u)  
	
	if (is.infinite(m))  
	{ 
		repeat
		{ if (L<=lower) break
			uni.slice.evals <<- uni.slice.evals + 1
			if (logv(L, n, lambda)<=logy) break
			L <- L - w
		}
		
		repeat
		{ if (R>=upper) break
			uni.slice.evals <<- uni.slice.evals + 1
			if (logv(R, n, lambda)<=logy) break
			R <- R + w
		}
	}
	
	else if (m>1)  
	{ 
		J <- floor(runif(1,0,m))
		K <- (m-1) - J
		
		while (J>0)
		{ if (L<=lower) break
			uni.slice.evals <<- uni.slice.evals + 1
			if (logv(L, n, lambda)<=logy) break
			L <- L - w
			J <- J - 1
		}
		
		while (K>0)
		{ if (R>=upper) break
			uni.slice.evals <<- uni.slice.evals + 1
			if (logv(R, n, lambda)<=logy) break
			R <- R + w
			K <- K - 1
		}
	}
	
	
	if (L<lower){
        L <- lower
	}
	if (R>upper){
        R <- upper
	}
	
	
	repeat{ 
		v1 <- runif(1,L,R)
		
		uni.slice.evals <<- uni.slice.evals + 1
		logv1 <- logv(v1, n, lambda)
		
		if (logv1>=logy) break
		
		if (v1>v0) 
		{ R <- v1
		}
		else 
		{ L <- v1
		}
	}
	
	return (v1)
	
}


mvt.mcmc <- function(X, niter, prior.lower.v, prior.upper.v,
                     prior.Mu0=rep(0, ncol(X)),
                     prior.Sigma0=diag(10000, ncol(X)),
                     prior.p=ncol(X), prior.V=diag(1, ncol(X)),
                     initial.v=NULL, initial.Sigma=NULL){

    if(!is.matrix(X) || ncol(X) < 2 || nrow(X) < 3){
        stop("The input observations of the 'mvt.mcmc' has to be a matrix
              of more than three rows and one column.")
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

	Sigma0.inv <- solve(prior.Sigma0)
	V.inv <- solve(prior.V)
	
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
	
	#assign matrice for saving results
	Mu.save <- matrix(0, niter, d)
	Sigma.save <- array(0, dim=c(d,d,niter))
	v.save <- rep(0, niter)
	
	for(j in 1:niter){
 		#update Mu
		A <- sum(lambda)*Sigma.inv + Sigma0.inv
		b <- Sigma.inv%*%colSums(X*lambda) + Sigma0.inv%*%prior.Mu0
		S <- solve(A)
		#browser()
		Mu <- rmvnorm(1, S%*%b, S)
	
		#update Sigma.inv
		S <- matrix(0, d, d)
		for(i in 1:n){
			S <- S + lambda[i]*t(X[i,]-Mu) %*% (X[i,]-Mu)
		}
   		Sigma.inv <- rwish(n+prior.p, solve(V.inv + S))
		
		#update lambda
		for(i in 1:n){
			shape <- (d+v) / 2
			rate <- v/2 + ((X[i,]-Mu) %*% Sigma.inv %*% t(X[i,]-Mu))/2
			lambda[i] <- rgamma(1, shape, rate)
		}
		
		#update v
		v <- v.slice(v, n, lambda, logv, w=1, m=Inf,
                     lower=prior.lower.v, upper=prior.upper.v, logv0=NULL)[1]
			
		#save results
		Mu.save[j,] <- Mu 
		Sigma.save[,,j] <- solve(Sigma.inv)
		v.save[j] <- v

	}
	list(Mu.save=Mu.save, Sigma.save=Sigma.save, v.save=v.save)
}


