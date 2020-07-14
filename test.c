#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include "nonparam.h"

int main (int argc, char **argv) {
  FILE *f;
  char buf[256];
  char *cp;
  float vals[1000];
  char groups[1000];
  int n = 0;
  int flag;
  int *sort_index;
  int *tie_counts;
  float *rank;

  float p;
  f = fopen(argv[1], "r");
  while (fgets(buf, 255, f)) {
    groups[n] = atoi(buf);
    cp = buf;
    while (!isspace(*cp)) cp++;
    while (isspace(*cp))cp++;
    vals[n] = atof(cp);
    n++;
  }

  sort_index = malloc(sizeof(int) * n);
  tie_counts = malloc(sizeof(int)*n);
  rank = malloc(sizeof(float)*n);
  if (sort_index == NULL || tie_counts == NULL || rank == NULL) {
    fprintf (stderr, "Could not alloc\n");
    exit(999);
  }

  p = nonparam_compar (vals, groups, n,
		       2, /* num_groups */
		       sort_index, rank, tie_counts, &flag);

  printf ("%f\t%d\n", p, flag);
}

