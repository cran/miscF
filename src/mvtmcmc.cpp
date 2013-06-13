#include <RcppArmadillo.h>

RcppExport SEXP logv(SEXP v_, SEXP n_, SEXP lambda_)
{    
    using namespace Rcpp;

	double v=as<double>(v_); 
    int n=as<int>(n_);
    NumericVector lambda(lambda_);

    double vh = v/2;
    double loglike=n*vh*log(vh) - n*lgamma(vh) + (vh-1)*sum(log(lambda)) - vh*sum(lambda);
	return(wrap(loglike));
		
}

RcppExport SEXP vslice(SEXP v0_, SEXP n_, SEXP lambda_, SEXP w_, SEXP m_, SEXP lower_, SEXP upper_){
    
	using namespace Rcpp;
    
    double v0=as<double>(v0_); 
    int n=as<int>(n_);
    NumericVector lambda(lambda_);
	double w=as<double>(w_);
    double m=as<double>(m_);
    double lower=as<double>(lower_);
    double upper=as<double>(upper_);
    
	double logv0 = as<double>(logv(wrap(v0), wrap(n), wrap(lambda)));
	double logy = logv0 - rexp(1)[0];
	double u = runif(1,0,w)[0];
	double L = v0 - u;
	double R = v0 + (w-u);
	
    int flag;
    double temp;
    double v1;
    
	if (!std::isfinite(m)){ 
		flag = 0;
        while (flag==0){
            if (L <= lower) flag=1;
            temp = as<double>(logv(wrap(L), wrap(n), wrap(lambda)));
            if (temp <= logy) flag=1;
		    L = L - w;
        }   
        
		
		flag = 0;
        while (flag==0){
            if (R >= upper) flag=1;
            temp = as<double>(logv(wrap(R), wrap(n), wrap(lambda)));
            if (temp <= logy) flag=1;
		    R = R + w;
        }
	}
    
	if (L < lower){
        L = lower;
	}
	if (R > upper){
        R = upper;
	}
	
	flag = 0;
	while(flag==0){ 
        double rrr = runif(1,L,R)[0];
 		v1 = rrr;
    	double logv1 = as<double>(logv(wrap(v1), wrap(n), wrap(lambda)));
		
		if (logv1 >= logy) flag=1;
		
		if (v1 > v0) 
		{ R = v1;
		}
		else 
		{ L = v1;
		}
	}

	return (wrap(v1));
	
}


RcppExport SEXP mvtmcmc(SEXP niter_, SEXP lambda_, SEXP SigmaInv_, SEXP Sigma0Inv_, SEXP X_, SEXP priorMu0_, SEXP priorP_, SEXP VInv_, SEXP v_, SEXP priorLowerV_, SEXP priorUpperV_){
 
	using namespace Rcpp;

    int niter = as<int>(niter_);
    arma::colvec lambda = as<arma::colvec>(lambda_);
	arma::mat SigmaInv = as<arma::mat>(SigmaInv_);
	arma::mat Sigma0Inv = as<arma::mat>(Sigma0Inv_);
	NumericMatrix X(X_);
    arma::colvec priorMu0 = as<arma::colvec>(priorMu0_);
    int priorP = as<int>(priorP_);
	arma::mat VInv = as<arma::mat>(VInv_);
    double v = as<double>(v_);
    double priorLowerV = as<double>(priorLowerV_);
    double priorUpperV = as<double>(priorUpperV_);

  
    int nrow = X.nrow(), ncol = X.ncol();
    int j, i, k;

    NumericVector MuSave(niter*ncol);
    NumericVector SigmaSave(niter*ncol*ncol);
    NumericVector vSave(niter);


	for(j=0; j <niter; j++){
		arma::mat A;
		arma::mat b;
		NumericVector Mu(ncol);

		// update Mu
 		A = sum(lambda)*SigmaInv + Sigma0Inv;
        NumericMatrix Xl(nrow, ncol);
        for (i=0; i<nrow; i++){
            for(k=0; k< ncol; k++){
                Xl(i,k) = lambda[i]*X(i,k);
            }
        }
        NumericVector colSum(ncol);
        for (i=0; i<ncol; i++){
            colSum[i] = sum(Xl(_,i));
        }
        arma::colvec colSumArma(colSum.begin(), colSum.size());
        b = SigmaInv*colSumArma + Sigma0Inv*priorMu0;

 		arma::mat S;
        S = arma::inv(A);
        Environment mvtnorm("package:mvtnorm");
        Function rmvnorm = mvtnorm["rmvnorm"];
        arma::mat mm = S*b;

        Mu = rmvnorm(1, mm, S);

        // update Sigma.inv
        arma::mat SS=arma::zeros<arma::mat>(ncol,ncol) ;
        for(i = 0; i < nrow; i++){
            NumericVector dd = sqrt(lambda[i])*(X(i,_)-Mu);
            arma::colvec ddArma(dd.begin(), dd.size());
            SS = SS + ddArma*arma::trans(ddArma);
        }

        Environment MCMCpack("package:MCMCpack");
        Function rwish = MCMCpack["rwish"];
        NumericMatrix Sig(ncol, ncol);
        Sig = rwish(nrow+priorP, arma::inv(VInv + SS));
	    SigmaInv = as<arma::mat>(Sig);

		// update lambda
        for(i = 0; i < nrow; i++){
            double shape = (ncol+v) / 2;
            NumericVector dd = X(i,_)-Mu;
            arma::colvec ddArma(dd.begin(), dd.size());
            double qq = arma::as_scalar(arma::trans(ddArma) * SigmaInv * ddArma);
            double rate = v/2 + qq/2;
            lambda[i] = R::rgamma(shape, 1/rate);
        }

        // update v
		double w=1;
		double m=INFINITY;
        v = as<double>(vslice(wrap(v), wrap(nrow), wrap(lambda), wrap(w), wrap(m), wrap(priorLowerV), wrap(priorUpperV)));

		// save results
		vSave(j) = v;
        for(i=0; i < ncol; i++){
			MuSave(i+j*ncol) = Mu[i];
		}
		arma::mat Sigma = arma::inv(SigmaInv);
        for(i=0; i < ncol; i++){
			for(k=0; k < ncol; k++){
				SigmaSave(i+k*ncol + j*ncol*ncol) = Sigma(i,k);
			}
		}

	}

	return List::create(Named( "Mu.save", MuSave),
                        Named( "Sigma.save", SigmaSave),
                        Named( "v.save", vSave));
}
