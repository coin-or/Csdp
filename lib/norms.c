/*
 * Calls the blas routine to compute a norm of a vector.
 */

#include <math.h>
#include "declarations.h"

double norm2(n,x)
     int n;
     double *x;
{
  double nrm;
  int incx=1;

  nrm=COIN_LAPACK_FUNC(dnrm2,DNRM2)(&n,x,&incx);
  
  return(nrm);
}

double norm1(n,x)
     int n;
     double *x;
{
  double nrm;
  int incx=1;

  nrm=COIN_LAPACK_FUNC(dasum,DASUM)(&n,x,&incx);
  
  return(nrm);
}

double norminf(n,x)
     int n;
     double *x;
{
  int i;
  double nrm;
  int incx=1;

  i=COIN_LAPACK_FUNC(idamax,IDAMAX)(&n,x,&incx);
  nrm=fabs(x[i-1]);
  
  return(nrm);
}



