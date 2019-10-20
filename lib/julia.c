// API used by the Julia wrapper at https://github.com/JuliaOpt/CSDP.jl

double getindex(block,i,j)
     struct blockrec block;
     int i;
     int j;
{
  if (i < 1 || i > block.blocksize)
    {
      printf("Invalid row index %d, it should be between 1 and %d\n", i, block.blocksize);
      exit(1);
    }
  if (j < 1 || j > block.blocksize)
    {
      printf("Invalid row index %d, it should be between 1 and %d\n", j, block.blocksize);
      exit(1);
    }
  switch (block.blockcategory)
	{
	  case DIAG:
        if (i == j)
          return block.data.vec[i];
        else
          return 0.0;
        break;
	  case MATRIX:
        return block.data.mat[ijtok(i,j,block.blocksize)];
        break;
	  default:
	    printf("getindex illegal block type %d\n", block.blockcategory);
	    exit(206);
	};

}

struct blockrec getblockrec(blockmat, blk)
     struct blockmatrix blockmat;
     int blk;
{
  if (blk < 1 || blk > blockmat.nblocks)
    {
      printf("Invalid block %d, it should be between 1 and %d\n", blk, blockmat.nblocks);
      exit(1);
    }
  return blockmat.blocks[blk];
}

// Structure used to look problem with a O(1) lookup from block number
// too the pointer of the corresponding block.
struct LoadingProblem
{
    int total_dimension;
    int num_constraints;
    struct blockmatrix *pC;
    double *a;
    struct sparseblock **constraint_block_lookup;
    struct constraintmatrix *constraints;
};

struct LoadingProblem* allocate_loading_prob(pC, block_dims, num_constraints, num_entries, printlevel)
    struct blockmatrix *pC;
    int *block_dims;
    int num_constraints;
    int *num_entries;
    int printlevel;
{
    int mat;
    int blk;
    int numentries;
    int length;
    struct sparseblock *sparse_block;

    if (pC->blocks < 0)
      {
        if (printlevel >= 1)
          printf("Invalid number of blocks in matrix C: %d\n", pC->blocks);
        exit(1);
      }
    pC->blocks = (struct blockrec *) safe_malloc((pC->nblocks + 1) * sizeof(struct blockrec));
    for (blk = pC->nblocks; blk > 0; blk--)
      {
	    pC->blocks[blk].blocksize = abs(block_dims[blk]);
        if (block_dims[blk] < 0)
          {
	        pC->blocks[blk].blockcategory=DIAG;
            length = 1 + abs(block_dims[blk]);
            pC->blocks[blk].data.vec = (double *) safe_malloc(length * sizeof(double));
          }
        else
          {
	        pC->blocks[blk].blockcategory=MATRIX;
            length = block_dims[blk] * block_dims[blk];
            pC->blocks[blk].data.mat = (double *) safe_malloc(length * sizeof(double));
          }
      }
    zero_mat(*pC);

    struct LoadingProblem *problem = (struct LoadingProblem *) safe_malloc(sizeof(struct LoadingProblem));

    problem->total_dimension = 0;
    for (blk = 1; blk <= pC->nblocks; blk++)
      problem->total_dimension += pC->blocks[blk].blocksize;

    if (num_constraints < 0)
      {
        if (printlevel >= 1)
          printf("Invalid number of constraints: %d\n", num_constraints);
        exit(1);
      }
    problem->num_constraints = num_constraints;

    problem->pC = pC;
    problem->a = (double *) safe_malloc((num_constraints + 1) * sizeof(double));

    problem->constraint_block_lookup = (struct sparseblock **) safe_malloc((num_constraints) * (pC->nblocks) * sizeof(struct sparseblock *));
    problem->constraints = (struct constraintmatrix *) safe_malloc((num_constraints + 1) * sizeof(struct constraintmatrix));

    for (mat = 1; mat <= num_constraints; mat++)
      {
        problem->constraints[mat].blocks = NULL;
        for (blk = pC->nblocks; blk > 0; blk--)
          {
            numentries = num_entries[ijtok(mat,blk,num_constraints)];
            if (numentries < 0)
              {
                if (printlevel >= 1)
                  printf("Invalid number of entries for constraint %d and block %d: %d\n", mat, blk, numentries);
                exit(1);
              }
            else if (numentries == 0)
              {
                problem->constraint_block_lookup[ijtok(mat,blk,num_constraints)] = NULL;
              }
            else
              {
                sparse_block = (struct sparseblock *) safe_malloc(sizeof(struct sparseblock));
                sparse_block->next = problem->constraints[mat].blocks;
                sparse_block->nextbyblock = NULL;
                sparse_block->entries = (double *) safe_malloc((numentries + 1) * sizeof(double));
                sparse_block->iindices = (int *) safe_malloc((numentries + 1) * sizeof(int));
                sparse_block->jindices = (int *) safe_malloc((numentries + 1) * sizeof(int));
                sparse_block->numentries = 0;
                sparse_block->blocknum = blk;
                sparse_block->blocksize = abs(block_dims[blk]);
                sparse_block->constraintnum = mat;
                sparse_block->issparse = 1;

                problem->constraint_block_lookup[ijtok(mat,blk,num_constraints)] = sparse_block;
                problem->constraints[mat].blocks = sparse_block;
              }
          }
      }

    return problem;
}

void free_loaded_prob(prob,X,y,Z)
     struct LoadingProblem* prob;
     struct blockmatrix X;
     double *y;
     struct blockmatrix Z;
{
  free_prob(prob->total_dimension,prob->num_constraints,*(prob->pC),prob->a,
            prob->constraints,X,y,Z);
}

void free_loading_prob(loading_prob)
    struct LoadingProblem* loading_prob;
{
    int mat;
    free(loading_prob->constraint_block_lookup);
    free(loading_prob);
}


void setconstant(loading_prob,mat,ent)
  struct LoadingProblem *loading_prob;
  int mat;
  double ent;
{
  loading_prob->a[mat] = ent;
}

int addentry(loading_prob,mat,blk,indexi,indexj,ent,allow_duplicates)
  struct LoadingProblem *loading_prob;
  int mat;
  int blk;
  int indexi;
  int indexj;
  double ent;
  int allow_duplicates;
{
  struct blockrec *blocks;
  struct sparseblock *sparse_block;
  int blksz;
  int itemp;
  int index;
  double *vec;

  if (ent == 0.0)
    return(0);

  /*
   * Arrange things so that indexi <= indexj.
   */

  if (indexi > indexj)
    {
      itemp=indexi;
      indexi=indexj;
      indexj=itemp;
    };

  if (mat == 0)
    {
      blocks = loading_prob->pC->blocks;
	  blksz = blocks[blk].blocksize;
	  if (blocks[blk].blockcategory == DIAG)
        {
          index = indexi;
          vec = blocks[blk].data.vec;
        }
      else
        {
          index = ijtok(indexi,indexj,blksz);
          vec = blocks[blk].data.mat;
        }
      if (!allow_duplicates && vec[index] != 0.0)
        return(1);
      vec[index] += ent;
	  if (indexi != indexj && blocks[blk].blockcategory == MATRIX)
        vec[ijtok(indexj,indexi,blksz)] += ent;
    }
  else
    {
      sparse_block = loading_prob->constraint_block_lookup[ijtok(mat,blk,loading_prob->num_constraints)];
      sparse_block->numentries += 1;
	  sparse_block->entries[sparse_block->numentries] = ent;
	  sparse_block->iindices[sparse_block->numentries] = indexi;
	  sparse_block->jindices[sparse_block->numentries] = indexj;
    }

  return(0);
}

void loaded_initsoln(problem,pX0,py0,pZ0)
     struct LoadingProblem* problem;
     struct blockmatrix *pX0;
     double **py0;
     struct blockmatrix *pZ0;
{
    initsoln(problem->total_dimension,problem->num_constraints,
             *problem->pC,problem->a,problem->constraints,pX0,py0,pZ0);
}
int loaded_sdp(problem,constant_offset,pX,py,pZ,ppobj,pdobj,printlevel,params)
     struct LoadingProblem* problem;
     double constant_offset;
     struct blockmatrix *pX;
     double **py;
     struct blockmatrix *pZ;
     double *ppobj;
     double *pdobj;
     int printlevel;
     struct paramstruc params;
{
    return parametrized_sdp(problem->total_dimension,problem->num_constraints,*problem->pC,
                            problem->a,problem->constraints,constant_offset,pX,py,pZ,ppobj,pdobj,printlevel,params);
}
