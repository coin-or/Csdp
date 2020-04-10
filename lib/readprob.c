/*
  Read in out a problem in SDPA sparse format.  Return 0 if ok, 1 if
  failure.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h>
#include "declarations.h"

void* safe_malloc();

#include "julia.c"

int load_prob_from_file();
int read_prob_size();
int load_prob();

FILE* sdpa_fopen();
void skip_to_end_of_line();
int safe_get_line();
int get_line();

int read_prob(fname,pn,pk,pC,pa,pconstraints,printlevel)
     char *fname;
     int *pn;
     int *pk;
     struct blockmatrix *pC;
     double **pa;
     struct constraintmatrix **pconstraints;
     int printlevel;

{
  int ret;
  int blk;
  struct LoadingProblem *loading_prob;

  ret = load_prob_from_file(fname,pC,&loading_prob,printlevel);
  if (ret != 0)
    return ret;

  *pn = loading_prob->total_dimension;
  *pk = loading_prob->num_constraints;
  *pa = loading_prob->a;
  *pconstraints = loading_prob->constraints;

  free_loading_prob(loading_prob);

  /*
   * If the printlevel is high, print out info on constraints and block
   * matrix structure.
   */
  if (printlevel >= 3)
    {
      printf("Block matrix structure.\n");
      for (blk=1; blk<=pC->nblocks; blk++)
	    {
	      if (pC->blocks[blk].blockcategory == DIAG)
	        printf("Block %d, DIAG, %d \n",blk,pC->blocks[blk].blocksize);
	      if (pC->blocks[blk].blockcategory == MATRIX)
	        printf("Block %d, MATRIX, %d \n",blk,pC->blocks[blk].blocksize);
	    };
    };

  return(0);
}

int load_prob_from_file(fname,pC,problem,printlevel)
     char *fname;
     struct blockmatrix *pC;
     struct LoadingProblem **problem;
     int printlevel;

{
  int ret;
  int buflen;
  char *buf;

  /*
   * This constant allows for up to 500000*20=10,000,000 byte long lines.
   * You might see a line this long in a problem with 500000 constraints.
   */

  buflen = 8000000;

  buf = (char *) safe_malloc(buflen*sizeof(char));

  /*
   * In the first pass through the file, determine the size parameters,
   * and allocate space for everything.
   */

  ret = read_prob_size(fname,pC,buf,buflen,problem,printlevel);
  if (ret != 0)
    {
      free(buf);
      return ret;
    }

  /*
   *  In the final pass through the file, fill in the actual data.
   */

  ret = load_prob(fname,buf,buflen,*problem,printlevel);
  free(buf);
  if (ret != 0)
    {
      return ret;
      free_loading_prob(*problem);
    }

  return(0);
}

int load_prob(fname,buf,buflen,loading_prob,printlevel)
     char *fname;
     char *buf;
     int buflen;
     struct LoadingProblem *loading_prob;
     int printlevel;
{
  char *ptr1;
  char *ptr2;
  FILE *fid;
  int ret;
  int mat;
  int blk;
  int indexi;
  int indexj;
  double ent;

  fid = sdpa_fopen(fname, printlevel);
  /*
   * Skip the number of constraints (primal variables in SDPA terminology).
   */
  ret = safe_get_line(fid,buf,buflen,"mDIM",printlevel);
  if (ret != 0) return(1);
  /*
   * Skip the number of blocks.
   */
  ret = safe_get_line(fid,buf,buflen,"nBLOCKS",printlevel);
  if (ret != 0) return(1);
  /*
   * Skip the block structure.
   */
  ret = safe_get_line(fid,buf,buflen,"block sizes",printlevel);
  if (ret != 0) return(1);

  /*
   *  Read in the right hand side values.
   */
  ret = safe_get_line(fid,buf,buflen,"values",printlevel);
  if (ret != 0) return(1);
  /*
   * Decode k numbers out of the buffer.  Put the results in
   * a.
   */
  ptr1 = buf;
  for (mat=1; mat<=loading_prob->num_constraints; mat++)
    {
      setconstant(loading_prob,mat,strtod(ptr1,&ptr2));
      /*
       * Check for a case where ptr2 didn't advance.  This indicates
       * a strtod failure.
       */
      if (ptr1 == ptr2)
        {
          if (printlevel >= 1)
            printf("Incorrect SDPA file. Can't read RHS values.\n");
          fclose(fid);
          return(1);
        };
      ptr1=ptr2;
    };

  /*
   * Now, read the actual entries.
   */
  ret=fscanf(fid,"%d %d %d %d %le ",&mat,&blk,&indexi,&indexj,&ent);
  do {
    /*
     * No need for sanity checking the second time around.
     */

    ret=addentry(loading_prob,mat,blk,indexi,indexj,ent,0);
    if (ret != 0)
      {
        if (printlevel >= 1)
          {
            printf("Incorrect SDPA file. Duplicate entry.\n");
            printf("mat=%d\n",mat);
            printf("blk=%d\n",blk);
            printf("indexi=%d\n",indexi);
            printf("indexj=%d\n",indexj);
          };

        fclose(fid);
        return(1);
      };
    ret=fscanf(fid,"%d %d %d %d %le ",&mat,&blk,&indexi,&indexj,&ent);
  } while (ret == 5);

  if ((ret != EOF) && (ret != 0))
    {
      if (printlevel >= 1)
        printf("Incorrect SDPA file. \n");
      fclose(fid);
      return(1);
    };


  fclose(fid);
  return(0);
}

int read_prob_size(fname,pC,buf,buflen,ploading_prob,printlevel)
     char *fname;
     struct blockmatrix *pC;
     char *buf;
     int buflen;
     struct LoadingProblem **ploading_prob;
     int printlevel;
{
  int ret;
  int blk;
  int mat;
  int indexi;
  int indexj;
  double ent;
  char *ptr1;
  char *ptr2;
  FILE *fid;
  int *isdiag;
  int *block_dims;
  int *num_entries;
  int num_constraints;

  fid = sdpa_fopen(fname, printlevel);

  /*
   * Get the number of constraints (primal variables in SDPA terminology)
   */
  ret = safe_get_line(fid,buf,buflen,"mDIM",printlevel);
  if (ret != 0) return(1);
  ret = sscanf(buf,"%d",&num_constraints);
  if ((ret != 1) || (num_constraints <= 0))
  {
    if (printlevel >= 1)
      printf("Incorrect SDPA file.  Couldn't read mDIM\n");
    fclose(fid);
    return(1);
  };

#ifndef BIT64
  /*
   * If operating in 32 bit mode, make sure that the dimension mDIM isn't
   * too big for 32 bits.  If we don't do this check, then integer overflow
   * won't be detected, and we'll allocate a bogus amount of storage for
   * O.
   */

  if (num_constraints > 23169)
    {
      if (printlevel >= 1)
        printf("This problem is too large to be solved in 32 bit mode!\n");
      exit(206);
    };
#endif

  /*
   * Read in the number of blocks.
   */
  ret = safe_get_line(fid,buf,buflen,"nBLOCKS",printlevel);
  if (ret != 0) return(1);
  ret = sscanf(buf,"%d",&(pC->nblocks));
  if ((ret != 1) || (pC->nblocks <= 0))
  {
    if (printlevel >= 1)
      printf("Incorrect SDPA file. Couldn't read nBLOCKS. \n");
    fclose(fid);
    return(1);
  }

  /*
   * And read the block structure.
   */
  block_dims = (int *) safe_malloc((pC->nblocks+1)*sizeof(int),printlevel);
  ret = safe_get_line(fid,buf,buflen,"block sizes",printlevel);
  if (ret != 0) return(1);
  /*
   * Decode nblocks numbers out of the buffer.  Put the results in
   * block_structure.
   */
  ptr1=buf;
  for (blk=1; blk<=pC->nblocks; blk++)
    {
      block_dims[blk] = strtol(ptr1,&ptr2,10);
      ptr1=ptr2;
    }

  /*
   *  Skip the right hand side values, it will be taken in to account in `load_prob`.
   */
  ret = safe_get_line(fid,buf,buflen,"values",printlevel);
  if (ret != 0) return(1);

  /*
   * Keep track of which blocks have off diagonal entries.
   */
  isdiag = (int *) safe_malloc((pC->nblocks+1)*sizeof(int),printlevel);
  for (blk=1; blk<=pC->nblocks; blk++)
    isdiag[blk] = 1;

  /*
   *  Now, loop through the entries,
   *  counting entries in the constraint matrices block by block.
   */

  num_entries = (int *) safe_malloc((num_constraints) * pC->nblocks * sizeof(int));
  for (mat = 1; mat <= num_constraints; mat++)
    for (blk = 1; blk <= pC->nblocks; blk++)
      num_entries[ijtok(mat,blk,num_constraints)] = 0;

  ret=fscanf(fid,"%d %d %d %d %le ",&mat,&blk,&indexi,&indexj,&ent);

  if (ret != 5)
    {
      if (printlevel >= 1)
        printf("Incorrect SDPA file. Return code from fscanf is %d, should be 5\n",ret);
      fclose(fid);
      free(isdiag);
      return(1);
    };

  do {

    /*
     * Check the validity of these values.
     */

    if ((mat < 0) || (mat > num_constraints) ||
	    (blk < 1) || (blk > pC->nblocks) ||
	    (indexi < 1) || (indexi > abs(block_dims[blk])) ||
	    (indexj < 1) || (indexj > abs(block_dims[blk])))
      {
        if (printlevel >= 1)
          printf("Incorrect SDPA file. Bad values in line: %d %d %d %d %e \n",
                 mat,blk,indexi,indexj,ent);
	    fclose(fid);
	    free(isdiag);
	    return(1);
      };

    /*
     * Mark this block as not diagonal if indexi!=indexj.
     */
    if (block_dims[blk] > 0 && indexi != indexj && ent != 0.0)
      isdiag[blk]=0;

    if (mat != 0)
      {
	    if (ent != 0.0)
	      num_entries[ijtok(mat,blk,num_constraints)] += 1;
      }
    else
      {
	    /*
	     * An entry in C. ignore it for now.
	     */
      };
    ret=fscanf(fid,"%d %d %d %d %le",&mat,&blk,&indexi,&indexj,&ent);
  } while (ret == 5);

  if ((ret != EOF) && (ret != 0))
    {
      if (printlevel >= 1)
        printf("Incorrect SDPA file, while reading entries.  ret=%d \n",ret);
      fclose(fid);
      free(isdiag);
      return(1);
    };

  /*
   * At this point, we'll stop to recognize whether any of the blocks
   * are "hidden LP blocks"  and correct the block type if needed.
   */

  for (blk=1; blk<=pC->nblocks; blk++)
    {
      if (block_dims[blk] > 1 && isdiag[blk] == 1)
	    {
	      /*
	       * We have a hidden diagonal block!
	       */
	      if (printlevel >= 2)
	        printf("Block %d is actually diagonal.\n",blk);
          block_dims[blk] = -block_dims[blk];
	    };
    };
  free(isdiag);

  *ploading_prob = allocate_loading_prob(pC, block_dims, num_constraints, num_entries, printlevel);

  free(block_dims);
  free(num_entries);

  fclose(fid);
  return(0);
}

void* safe_malloc(size, printlevel)
    int size;
    int printlevel;
{
  void *ptr = malloc(size);
  if (ptr == NULL)
    {
      if (printlevel >= 1)
        printf("Storage allocation of %d bytes failed!\n", size);
      exit(205);
    }
  return ptr;
}

FILE* sdpa_fopen(fname, printlevel)
     char *fname;
     int printlevel;
{
  char c;
  FILE *fid;
  fid=fopen(fname,"r");

  if (fid == (FILE *) NULL)
    {
      if (printlevel >= 1)
        printf("Couldn't open problem file for reading! \n");
      exit(201);
    };

  /*
   * First, read through the comment lines.
   */

  c=getc(fid);
  while ((c == '"') || (c == '*'))
    {
      skip_to_end_of_line(fid);
      c=getc(fid);
    };

  ungetc(c,fid);

  return fid;
}

/*
 *  This routine skips to the end of the current line of input from the
 *  file fid.
 */

void skip_to_end_of_line(fid)
     FILE *fid;
{
  char c;

  c=getc(fid);
  while (c != '\n')
    c=getc(fid);
}

int safe_get_line(fid,buffer,bufsiz,name,printlevel)
     FILE *fid;
     char *buffer;
     int bufsiz;
     char *name;
     int printlevel;
{
  int ret;
  ret = get_line(fid,buffer,bufsiz);
  if (ret != 0)
    {
      if (printlevel >= 1)
        printf("Incorrect SDPA file. Can't read %s.\n", name);
      fclose(fid);
    }
  return ret;
}

/*
 *  This routine reads a line of input into a buffer, and translates all
 *  occurences of "," "{" "}" "(" ")" to spaces.
 *
 */

int get_line(fid,buffer,bufsiz)
     FILE *fid;
     char *buffer;
     int bufsiz;
{
  int k;
  char c;

  if (fgets(buffer,bufsiz-1,fid) != NULL)
    {
      c=buffer[0];
      k=0;
      while ((c != '\n')  && (c!=0) && (k<bufsiz))
        {
          if (c==',')
            buffer[k]=' ';
          if (c=='(')
            buffer[k]=' ';
          if (c==')')
            buffer[k]=' ';
          if (c=='{')
            buffer[k]=' ';
          if (c=='}')
            buffer[k]=' ';
          k=k+1;
          c=buffer[k];
        };

      /*
       * return an error if the line is longer than the buffer!
       */

      if (k<bufsiz-5)
        return(0);
      else
        return(1);
    }
  else
    return(2);
}
