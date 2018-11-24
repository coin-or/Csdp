/*
  Compute the dual objective function value dobj=a'y.  
  */

#include "declarations.h"

double calc_dobj(k,a,y,constant_offset)
     int k;
     double *a;
     double *y;
     double constant_offset;
{
  double s;
  int incx=1;

  s=0.0;

  s=s+COIN_LAPACK_FUNC(ddot,DDOT)(&k,a+1,&incx,y+1,&incx);

  return(s+constant_offset);
  
}

