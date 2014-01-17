#include "Rmath.h"
#include "R.h"

double soft( double x, double t){
  if(fabs(x) - t < 0)
    return 0;
  else if(x > 0)
    return x - t;
  else
    return x + t;
}

void lasso(double *r, double *X, double *y, double *beta, int *p, int *n, double *lambda, int *maxit, double *eps){
  double b, betanew, maxdiff;
  int its = 0;
  do{
    maxdiff = 0;
    for(int i=0; i < *p; i++){  
      b=0;
      for(int k=0; k < *n; k++){
        r[k] += beta[i]*X[k + i*(*n)];
        b += r[k]*X[k + i*(*n)];
      }
      b=b/(*n);
      //Rprintf("\niterations %i:variable \t%i: beta (pre-thresh)\t%e: maxdiff \t%e",its, i, b, maxdiff);
      betanew = soft(b, *lambda);
      maxdiff = fmax(fabs(betanew - beta[i]),maxdiff);
      beta[i] = betanew;
      //Rprintf("\niterations %i:variable \t%i: beta (post-thresh)\t%e: maxdiff \t%e",its, i, betanew, maxdiff);
      for(int k=0; k < *n; k++){
        r[k] -= beta[i]*X[k + i*(*n)];
      }
    }
    its++;
  } while(maxdiff > *eps & its < *maxit);
  return;
}




void Wlasso(double *w, double *r, double *X, double *y, double *beta, int *p, int *n, double *lambda, int *maxit, double *eps){
  double b, betanew, maxdiff, xwx;
  int its = 0;
  do{
    maxdiff = 0;
    for(int i=0; i < *p; i++){  
      xwx = 0;
      b=0;
      for(int k=0; k < *n; k++){
        r[k] += beta[i]*X[k + i*(*n)];
        b += r[k]*w[k]*X[k + i*(*n)];
	xwx += X[k + i*(*n)]*w[k]*X[k + i*(*n)];
      }
      betanew = soft(b, lambda[i])/xwx;
      maxdiff = fmax(fabs(betanew - beta[i]),maxdiff);
      beta[i] = betanew;
      for(int k=0; k < *n; k++){
        r[k] -= beta[i]*X[k + i*(*n)];
      }
    }
    its++;
  } while(maxdiff > *eps & its < *maxit);
  return;
}

//Rprintf("\n%i:\t%i:\t%e:\t%e",its, i, b, maxdiff);

