#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <math.h>

SEXP curvePolynomialGetDM(SEXP sl0, SEXP sl, SEXP sx, SEXP sxknots)
{
    int l0 = asInteger(sl0);
    int l = asInteger(sl);
    double *x = REAL(sx);
	double *xknots = REAL(sxknots);
    
	
	int I = LENGTH(sx);
	int nknot = LENGTH(sxknots);
	
	int d = l - l0;

	SEXP val;
	PROTECT(val = allocMatrix(REALSXP, I, l + (d+1)*nknot));
	double *dm = REAL(val);

	int *larger = (int *) R_alloc(I*nknot, sizeof(int));

	
	int i, j, k;
	
	for (i = 0; i < I; i++) {
		for (k = 0; k < nknot ; k++) {
			if (x[i] > xknots[k]) {
				larger[i + k*I] = 1;
			}
			else{
				larger[i + k*I] = 0;
			}
		}
	}

	
	for (i = 0; i < I; i++) {
		for(j = 0; j < l; j++) {
			dm[i + I * j] = pow(x[i], j+1);
		}
	}

	for (i = 0; i < I; i++) {
		for (k = 0; k < nknot ; k++) {
			for (j = l0; j < l+1; j++) {
				dm[i + I*(j+d + k*(d+1))] = pow(x[i]-xknots[k], j) * larger[i + k*I];
			}
		}
    }

  
    UNPROTECT(1);
    return val;
}

