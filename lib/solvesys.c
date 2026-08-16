/*
  Solve a system of equations using the Cholesky factorization of A.  
  Note that we assume that A is positive definite and that A has already
  been factored.
*/

#include <stdlib.h>
#include <stdio.h>
#include "declarations.h"

int solvesys(m,ldam,A,rhs)
     int m;
     int ldam;
     double *A;
     double *rhs;
{
  int incx;
  int info;

  incx=1;

  csdp_dpotrs("U",&m,&incx,A,&ldam,rhs+1,&ldam,&info);

	   if (info != 0)
	     {
	       return(6);
	     };

	   return(0);
}


