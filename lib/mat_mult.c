/*
 *  Compute C=scale1*A*B+scale2*C.  
 *  Note that C must consist of dense matrix and vector blocks- no sparse
 *  blocks or eye's or other special cases.
 *
 *  A and B can have blocks of all supported types.  Unsupported types
 *  generate exit(206).
 *
 *  It is assumed that all three matrices are of compatible block strucutre.
 *
 *  We use dgemm to do the actual work on normal dense blocks.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include "declarations.h"

void mat_mult(scale1,scale2,A,B,C)
     double scale1,scale2;
     struct blockmatrix A,B,C;
{
  int blk,i,n;
  double *ap;
  double *bp;
  double *cp;

  /*
   * In theory, the BLAS ensures that if scale2=0, then C will not be 
   * accessed before being written to.  Some poor implementations of the
   * BLAS would generate NaN results if C contained NaNs, and scale2=0.0.
   * This doesn't appear to be conforming to the reference implementation, 
   * and I'm not aware of any current BLAS implementations with this problem, 
   * so I'm commenting out this bit of code.  It's here only in case this
   * problem arises again.


         if (scale2 == 0.0)
            zero_mat(C);

   */

  /*
   * Work through the blocks one at a time.
   */
  
  for (blk=1; blk<=A.nblocks; blk++)
    {
      switch (A.blocks[blk].blockcategory) 
	{
	case DIAG:
          if (scale2 != 0.0)
	    {
	      for (i=1; i<=A.blocks[blk].blocksize; i++)
		{
		  C.blocks[blk].data.vec[i]=scale1*A.blocks[blk].data.vec[i]
		    *B.blocks[blk].data.vec[i]
		    +scale2*C.blocks[blk].data.vec[i];
		};
	    }
	  else
	    {
	      for (i=1; i<=A.blocks[blk].blocksize; i++)
		{
		  C.blocks[blk].data.vec[i]=scale1*A.blocks[blk].data.vec[i]
		    *B.blocks[blk].data.vec[i];
		};
	    };
	  break;
	case MATRIX:
	  /*
	   * Call dgemm to do the matrix multiplication.
	   */
	  
	  n=A.blocks[blk].blocksize;
	  ap=A.blocks[blk].data.mat;
	  bp=B.blocks[blk].data.mat;
	  cp=C.blocks[blk].data.mat;

	  mat_mult_raw(n,scale1,scale2,ap,bp,cp);
	  break;
	default:
	  printf("mat_mult illegal block type!\n");
	  exit(206);
	};

    };
}

void mat_mult_raw(n,scale1,scale2,ap,bp,cp)
     int n;
     double scale1;
     double scale2;
     double *ap;
     double *bp;
     double *cp;
{

#ifdef HIDDENSTRLEN
  dgemm_("N","N",&n,&n,&n,&scale1,ap,&n,bp,&n,&scale2,cp,&n,1,1);
#else
  CSDP_LAPACK_FUNC(dgemm,DGEMM)("N","N",&n,&n,&n,&scale1,ap,&n,bp,&n,&scale2,cp,&n);
#endif
}

#ifdef USEATLAS


void mat_mult_rawatlas(n,scale1,scale2,ap,bp,cp)
     int n;
     double scale1;
     double scale2;
     double *ap;
     double *bp;
     double *cp;
{
enum CBLAS_ORDER {CblasRowMajor=101, CblasColMajor=102 };
enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113,
		      AtlasConj=114};

void  ATL_dgemm(enum CBLAS_TRANSPOSE transa,enum CBLAS_TRANSPOSE transb,
                int m,int n, int k,double scale1,double *ap,int l,double *bp,
                int p,double scale2,double *cp,int q);

 ATL_dgemm(CblasNoTrans,CblasNoTrans,n,n,n,scale1,ap,n,bp,n,scale2,cp,n);
}

#endif

