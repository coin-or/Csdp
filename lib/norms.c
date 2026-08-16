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

  nrm=csdp_dnrm2(&n,x,&incx);
  
  return(nrm);
}

double norm1(n,x)
     int n;
     double *x;
{
  double nrm;
  int incx=1;

  nrm=csdp_dasum(&n,x,&incx);
  
  return(nrm);
}

double norminf(n,x)
     int n;
     double *x;
{
  int i;
  double nrm;
  int incx=1;

  i=csdp_idamax(&n,x,&incx);

  return(nrm);
}



