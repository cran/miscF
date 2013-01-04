#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>

SEXP spatialLmeUpdateMu(SEXP sI, SEXP sG, SEXP scV, SEXP ssV, SEXP smu0,
			  SEXP sbetaBar, SEXP salpha, SEXP ssigma2inv, SEXP slambda2inv)
{
	int I = asInteger(sI);
	int G = asInteger(sG);
    int *cV = INTEGER(scV);
    int sV = asInteger(ssV);
	double *mu0 = REAL(smu0);
    double *betaBar = REAL(sbetaBar);
	double *alpha = REAL(salpha);
	double *sigma2inv = REAL(ssigma2inv);
    double *lambda2inv = REAL(slambda2inv);

    SEXP val = allocVector(REALSXP, sV);
    double *mu = REAL(val);
    
	double *alphaBar = (double *) R_alloc(G, sizeof(double));
	
	int i, g, v;
	
	GetRNGstate();
		
    for (g = 0; g < G; g++) {
        alphaBar[g] = 0;
		for (i = 0; i < I; i++) {
			alphaBar[g] += alpha[g + i*G];
		}
		alphaBar[g] = alphaBar[g] / I;
    }
	
	for(g = 0; g < G; g++) {
		double m, var, s, isig;
		for( v = cV[g]; v < cV[g+1]; v++) {
			m = betaBar[v] - alphaBar[g];
			isig = I*sigma2inv[g];
			m = m * isig;
			m = m + (mu0[g]*lambda2inv[g]);
			var = 1 / (isig + lambda2inv[g]);
			m = m * var;
			s = sqrt(var);
			mu[v] = norm_rand() * s + m;
		}
				
	}

    PutRNGstate();

   return val;
}

