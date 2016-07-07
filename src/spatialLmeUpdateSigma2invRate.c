#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>

SEXP spatialLmeUpdateSigma2invRate(SEXP sspatialMat, SEXP sI, SEXP sG, SEXP scV, 
						 SEXP ssumV, SEXP smu, SEXP salpha)
{
	int I = asInteger(sI);
	int G = asInteger(sG);
    int *cV = INTEGER(scV);
    int sumV = asInteger(ssumV);
    double *spatialMat = REAL(sspatialMat);
	double *mu = REAL(smu);
	double *alpha = REAL(salpha);

    SEXP val = allocVector(REALSXP, G);
    double *rate = REAL(val);
    

	int i, g, v;
	
	double r;

	for (g = 0; g < G; g++) {
		rate[g] = 0;
		for (i = 0; i < I; i++) {
			for( v = cV[g]; v < cV[g+1]; v++) { 
				r = spatialMat[i*sumV + v] - mu[v] - alpha[i*G + g];
				rate[g] += r * r;
			}
		}
	}
	
    return val;
}

