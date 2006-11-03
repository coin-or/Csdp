
#include <stdlib.h>
#include <stdio.h>
#include "declarations.h"

struct entry {
#ifndef NOSHORTS
  unsigned short indexi;
  unsigned short indexj;
  double entry;
#else
  int indexi;
  int indexj;
  double entry;
#endif
};

int mycompare(p1,p2)
     struct entry *p1;
     struct entry *p2;
{
  if (p1->indexi < p2->indexi)
    {
      return(-1);
    };
  if (p1->indexi == p2->indexi)
    {
      if (p1->indexj < p2->indexj)
	{
	  return(-1);
	};
      if (p1->indexj == p2->indexj)
	{
	  return(0);
	};
      if (p1->indexj > p2->indexj)
	{
	  return(1);
	};
    };

  /*
   * If we get here, then p1->indexi > p2->indexi
   */
  
  return(1);

}

void sort_entries(k,C,constraints)
     int k;
     struct blockmatrix C;
     struct constraintmatrix *constraints;
{
  int i,j;
  struct sparseblock *ptr;
  int maxentries;
  struct entry *entries;

  /*
   * First, find out who has the most entries.
   */

  maxentries=0;
  for (i=1; i<=k; i++)
    {
      ptr=constraints[i].blocks;
      while (ptr != NULL)
	{
	  if (ptr->numentries > maxentries)
	    maxentries=ptr->numentries;
	  ptr=ptr->next;
	};
    };

  /*
   * Allocate space for entries.
   */

  entries=(struct entry *)malloc(maxentries*sizeof(struct entry));
  if (entries==NULL)
    {
      printf("Storage allocation failed in sortentries.\n");
      exit(10);
    };
  
  for (i=1; i<=k; i++)
    {
      ptr=constraints[i].blocks;
      while (ptr != NULL)
	{
	  /*
	   * Copy in
	   */

	  for (j=1; j<=ptr->numentries; j++)
	    {
	      entries[j-1].indexi=ptr->iindices[j];
	      entries[j-1].indexj=ptr->jindices[j];
	      entries[j-1].entry=ptr->entries[j];
	    };

	  /*
	   * Sort
	   */
	  
	  qsort(entries,(size_t)ptr->numentries,sizeof(struct entry),
		&mycompare);

	  /*
	   * Copy out.
	   */

	  for (j=1; j<=ptr->numentries; j++)
	    {
	      ptr->iindices[j]=entries[j-1].indexi;
	      ptr->jindices[j]=entries[j-1].indexj;
	      ptr->entries[j]=entries[j-1].entry;
	    };


	  ptr=ptr->next;
	};
    }; /* end i */

  free(entries);


}
