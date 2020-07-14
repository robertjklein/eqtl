/*
 * structs.h
 *
 * Main data structures for program
 *
 * Robert J. Klein
 * July 8, 2010
 */

#ifndef _structs_h
#define _structs_h

#define MAXP 0.05
#define ALPHA 0.05
#define FDR_ALPHA 0.1
#define THRESHOLD 0.05 

#ifndef VERSION
#define VERSION "0.0"
#endif

#ifndef COPYRIGHT
#define COPYRIGHT "Copyringt (C) 2010 Memorial Sloan-Kettering Cancer Center"
#endif

#ifndef LICENSE
#define LICENSE "Licensed under the terms of the GNU General Public License (GPL).  See file \nLICENSE for more information.\n"
#endif

typedef struct _snp_t {
  char chr;
  int pos;
  int rs;
  char *gt;
  char **id_list;
  int num_indivs;
  int num_snps;
  int num_groups;
  struct _snp_t *next;
} snp_t;

typedef struct _phen_t {
  char *name;
  float *values;
  char chr;
  int start;
  int stop;
  struct _phen_t *next;
} phen_t;

typedef struct _result_t {
  snp_t *snp;
  phen_t *phen;
  double p;
  int flag;
  char good_for_cis;
} result_t;

#endif
