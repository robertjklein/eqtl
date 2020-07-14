/*
 * nonparam.h
 *
 * Routines for simple nonparametric statistics for eQTL analysis.
 * Specifically, does Kruskal-Wallis for 3 groups or Mann-Whitney for
 * two groups.  Otherwise, fails.  Assumes needed array space already
 * allocated and passed in.
 *
 * Robert J. Klein
 * June 19, 2010
 */

#ifndef _nonparam_h
#define _nonparam_h

double nonparam_compar (float *vals, char *groups, int n, int num_groups\
		       , int *sort_index, float *rank, int *tie_counts, int *flag);

#endif
