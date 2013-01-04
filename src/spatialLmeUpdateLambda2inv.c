#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>

SEXP spatialLmeUpdateLambda2inv(SEXP sG, SEXP sV, SEXP scV, 
								SEXP smu, SEXP smu0, SEXP sc0, SEXP sd0)
{
	int G = asInteger(sG);
    int *V = INTEGER(sV);
	int *cV = INTEGER(scV);
	double *mu = REAL(smu);
	double *mu0 = REAL(smu0);
	double c0 = asReal(sc0);
	double d0 = asReal(sd0);

    SEXP val = allocVector(REALSXP, G);
    double *lambda2inv = REAL(val);
    

	int g, v;
	double rate, shape;

	GetRNGstate();
	
	for (g = 0; g < G; g++) {
		rate = 0;
		double m = mu0[g];
		for(v = cV[g]; v < cV[g+1]; v++) { 
			double r = mu[v] - m;
  			rate += r * r;
		}
		rate = d0 + 0.5*rate;
		shape = c0 + V[g]/2;
		lambda2inv[g] = rgamma(shape, 1/rate);
	}
	
    PutRNGstate();
	
    return val;
}

