/*
 * User exit routine for psd.  This version of the routine simply returns
 * 0, so that psd will not stop.
 */

#include "declarations.h"

int user_exit(n,k,C,a,dobj,pobj,constant_offset,constraints,X,y,Z,params)
     int n;
     int k;
     struct blockmatrix C;
     double *a;
     double dobj;
     double pobj;
     double constant_offset;
     struct constraintmatrix *constraints;
     struct blockmatrix X;
     double *y;
     struct blockmatrix Z;     
     struct paramstruc params;
{
  return(0);
}

