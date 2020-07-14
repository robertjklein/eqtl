/*
 * nonparam.h
 *
 * Routines for simple nonparametric statistics
 *
 * Robert J. Klein
 * June 19, 2010
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_cdf.h>
#include "nonparam.h"

/* 
 * Given a list of values and groupings, does kruskal-wallis
 * if 3 groups or Mann-Whitney if 2.  Returns p-value.  Includes
 * a small value flag if p-value is to be believed.  Computes rank sum for both
 * in main routine b/c code is the same.
 */
double nonparam_compar (float *vals, char *groups, int n, int num_groups, int *sort_index, float *rank, int *tie_counts, int *flag) {
  float H, U;
  int i, j, sum, tot_ties;
  float avg_rank;
  float rank_sum[3];
  float n_i[3];

  /* First sort the data */
  for (i=0; i<n; i++) {
    sort_index[i] = i;
  }
  int sort_func (const void *a, const void *b) {
    int i, j;

    i = *((int *)a);
    j = *((int *)b);

    if (vals[i] < vals[j]) {
      return(-1);
    } else if (vals[i] > vals[j]) {
      return(1);
    } else {
      return (0);
    }
  }
  qsort (sort_index, n, sizeof(int), &sort_func);

  /* Remove missing indivs (group = 127) from sort index */
  i=0; j=0; 
  while (i<n) {
    if ((int)groups[sort_index[i]] != 127) {
      sort_index[j] = sort_index[i];
      j++;
    }
    i++;
  }
  n = j;

  /* Set up tie counts */
  for (i=0; i<n; i++) {
    tie_counts[i] = 0;
  }
  tot_ties = 0;

  /* Compute ranks */
  for (i=0; i<n; i++) {
    if (vals[sort_index[i]] != vals[sort_index[i+1]]) {
      rank[sort_index[i]] = i + 1.;
    } else {
      sum = 0;
      j = i;
      while (j<n && vals[sort_index[i]] == vals[sort_index[j]]) {
	sum += (j+1);
	j++;
      }
      avg_rank = (float)sum / (j-i);
      tie_counts[tot_ties] = j-i;
      tot_ties++;
      while (i<j) {
	rank[sort_index[i]] = avg_rank;
	i++;
      }
      i--;
    }
  }

  /* Compute rank sums and n_i */
  for (i=0; i<3; i++) {
    n_i[i] = 0;
    rank_sum[i] = 0.;
  }
  for (i=0; i<n; i++) {
    n_i[(int)groups[sort_index[i]]]++;
    rank_sum[(int)groups[sort_index[i]]] += rank[sort_index[i]];
  }

  *flag = 0;

  /* Mann-Whitney for 2 groups */
  if (num_groups == 2) {
    H = rank_sum[0] - 0.5*(float)n_i[0]*(n_i[0]+1.);
    if (n_i[0]*n_i[1] - H > H) {
      U = n_i[0]*n_i[1] - H;
    } else {
      U = H;
    }
    /* subtract out mean */
    U -= (0.5*n_i[0]*n_i[1]);
    
    /* Use H for std. deviation here */
    if (tot_ties > 0) {
      for (i=0; i<tot_ties; i++) {
	H += ((float)(tie_counts[i]*(tie_counts[i]+1)*(tie_counts[i]-1)))/12.;
      }
      H *= ((float)(n_i[0]*n_i[1]))/((float)(n*(n-1)));
      H *= -1.;
    } else {
      H = 0.;
    }
    H += ((float)(n_i[0]*n_i[1]*(n+1)))/12.;
    U /= sqrtf(H);
    *flag = -1;

    return(2.*gsl_cdf_ugaussian_Q(fabs((double)U)));

    /* Kruskal-Wallis for 3 groups */
  } else if (num_groups == 3) {
    *flag = 1;
    H = (12 * 
      (((rank_sum[0]*rank_sum[0])/(float)n_i[0]) +
      ((rank_sum[1]*rank_sum[1])/(float)n_i[1]) +
       ((rank_sum[2]*rank_sum[2])/(float)n_i[2])) /
	 ((float)(n)*(float)(n+1))) - 
      3 *(n+1);
    if (tot_ties > 0) {
      sum = 0;
      for (i=0; i<tot_ties; i++) {
	sum += (tie_counts[i]*tie_counts[i]*tie_counts[i]-tie_counts[i]);
      }
      H = H/(1. - (float)sum/((float)(n*n*n-n)));
    }
    for (i=0; i<3; i++) {
      if (n_i[i] < 5) *flag = 2;
    }
    return(exp(-0.5*(double)H));
  } else {
    *flag = 2;
    return(-1.);
  }
}
