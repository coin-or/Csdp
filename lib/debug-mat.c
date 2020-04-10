#include <stdio.h>
#include "blockmat.h"


void printb(struct blockrec b) {
  int i;
  int size;
  printf("blockrec:\n");
  printf(" blockcategory: %d\n", (int) b.blockcategory);
  printf(" blocksize:     %d\n", b.blocksize);
  printf(" data.vec:      %p\n", b.data.vec);
  if (b.blocksize <= 16) {
    switch (b.blockcategory) {
    case MATRIX:
      size = b.blocksize * b.blocksize;
      for (i = 0; i < size; i++)
        printf("  %f\n", b.data.vec[i]);
      break;
    case DIAG:
      for (i = 1; i <= b.blocksize; i++)
        printf("  %f\n", b.data.vec[i]);
      break;
    default:
      fprintf(stderr, "UNKNOWN blockcat: %d",
              (int) b.blockcategory);
    }
  }
}


void printb_(struct blockrec *b) {
  printb(*b);
}


void print_sizeof() {
  printf("sizeof(void*)                 %d\n",
         (int) sizeof(void*));
  printf("sizeof(int)                   %d\n",
         (int) sizeof(int));
  printf("sizeof(enum blockcat)         %d\n",
         (int) sizeof(enum blockcat));
  printf("sizeof(struct blockrec)       %d\n",
         (int) sizeof(struct blockrec));
  printf("sizeof(struct blockmatrix)    %d\n",
         (int) sizeof(struct blockmatrix));
  printf("sizeof(struct sparseblock)    %d\n",
         (int) sizeof(struct sparseblock));
}



void printm(struct blockmatrix A) {
  int blk;
  printf("A.nblocks = %d\n", A.nblocks);
  printf("A.blocks  = %p\n", A.blocks);
  for (blk=1; blk <= A.nblocks; blk++) {
    printf("block[%d]: %p\n", blk, &A.blocks[blk]);
    printb(A.blocks[blk]);
  }

}

FILE *fid;

void print_sparse_block(struct sparseblock *b) {
  int i;
  fid = stdout;
  fprintf(fid, "\n* Printing block: %p\n", b);
  if (b == NULL)
    return;
  fprintf(fid, " next: %p\n", b->next);
  fprintf(fid, " nextbyblock: %p\n", b->nextbyblock);
  fprintf(fid, " constraintnum: %d\n", b->constraintnum);
  fprintf(fid, " blocknum: %d\n", b->blocknum);
  fprintf(fid, " blocksize: %d\n", b->blocksize);
  fprintf(fid, " numentries: %d\n", b->numentries);
  if (b->blocksize <= 30) {
    for (i = 1; i <= b->numentries; i++) {
      fprintf(fid, "  block[%d, %d] = %f\n",
              b->iindices[i],
              b->jindices[i],
              b->entries[i]);
    }
  } else {
    fprintf(fid, " TOO LARGE\n");
  }
  print_sparse_block(b->next);
}

void  print_constraints(int k,
                        struct constraintmatrix *constraints)
{
  int i, j;
  fid = stdout;
  struct sparseblock *p;

  fprintf(fid, "constraints == %p\n", constraints);

  for (i=1; i<=k; i++)
    {
      fprintf(fid, "\n\nprinting constraints[%d].blocks\n", i);
      p=constraints[i].blocks;
      while (p != NULL)
        {
          fprintf(fid, "p == %p\n\n", p);
          print_sparse_block(p);
          fprintf(fid, "\n");
          for (j=1; j<=p->numentries; j++)
            {
              fprintf(fid, "i=%d, j=%d\n", i, j);
              fprintf(fid,"%d %d %d %d %f \n",
                      i,
                      p->blocknum,
                      p->iindices[j],
                      p->jindices[j],
                      p->entries[j]);
            };
          p=p->next;
        };
    };
}
