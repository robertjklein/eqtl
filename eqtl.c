/*
 * eqtl.c
 *
 * Program for doing genome-wide eQTL analysis (all SNPs x all probes, in 
 * theory).
 *
 * Robert J. Klein
 * June 18, 2010
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "squid.h"
#include "sqfuncs.h"

#include "structs.h"
#include "eqtlio.h"
#include "nonparam.h"
#include "regress.h"

static char banner[] = "eqtl -- performs genome wide eQTL analysis\n";

static char usage[] = "\
Usage: eqtl [-options] <PLINK prefix> <gene list> <expression directory>\n\
  Available optiosn are:\n\
  -h      : help; print brief help on version and udage\n\
  -c      : Look for cis-eQTLs only\n\
";

static char experts[] = "\
   --test <s>     : Specifies test to do.  Options are kw [krusal-wallis] or reg [linear regression]\n\
   --qnorm        : Quantile normalize the expression data\n\
   --dist <kb>    : Kilobases to do cis search in\n\
";

static struct opt_s OPTIONS[] = {
  { "-h", TRUE, sqdARG_NONE },
  { "-c", TRUE, sqdARG_NONE },
  { "--test", FALSE, sqdARG_STRING },
  { "--qnorm", FALSE, sqdARG_NONE },
  { "--dist", FALSE, sqdARG_INT }
};
#define NOPTIONS (sizeof(OPTIONS)/sizeof(struct opt_s))

int check_cis (snp_t *snp, phen_t *phen, int maxdist) {
  int retval = 0;
  
  if (snp->chr == phen->chr &&
      (((snp->pos < phen->start) && (phen->start - snp->pos <= maxdist)) ||
       ((snp->pos > phen->stop) && (snp->pos - phen->stop <= maxdist)) ||
       ((snp->pos >= phen->start) && (snp->pos <= phen->stop)))) {
    retval = 1;
  }
  return(retval);
}

int result_sort_func (const void *a, const void *b) {    
  result_t **i;
  result_t **j;

  i = (result_t **)a;
  j = (result_t **)b;

  if ((*i)->p < (*j)->p) {
    return(-1);
  } else if ((*i)->p > (*j)->p) {
    return(1);
  } else {
    return (0);
  }
}

result_t **get_results (snp_t *genotypes, phen_t *phenotypes, long long *tot_results_r, long long *total_tests_r, long long *total_cis_tests_r, int test_type, int cis_only, int maxdist) { 
  long long estimated_results;
  int phen_count = 0;
  int iteration = 0;

  phen_t *cur_phen;
  snp_t *cur_snp;

  float p = -1.0;
  int flag;
  int is_cis;
  result_t **results;

  long long tot_results = 0;
  long long total_tests = 0;
  long long total_cis_tests = 0;

  int num_indivs;

  int *sort_index = NULL;
  float * rank = NULL;
  int *tie_counts = NULL;

  /* Prepare for analysis */
  num_indivs = genotypes->num_indivs;

  if (test_type == 0) {
    sort_index = MallocOrDie(sizeof(int)*num_indivs);
    tie_counts = MallocOrDie(sizeof(int)*num_indivs);
    rank = MallocOrDie(sizeof(float)*num_indivs);
  }

  /* Count phenotypes and let us know how many tests */
  for (cur_phen=phenotypes; cur_phen != NULL; cur_phen = cur_phen->next) {
    phen_count++;
  }
  printf ("There are %d snps in %d phenotypes tested in %d individuals\n", genotypes->num_snps, phen_count, num_indivs);

  /* Get estimated results memory */
  estimated_results = 4*((long long)(MAXP * genotypes->num_snps * phen_count));
  results = malloc(sizeof(result_t *)*estimated_results);
  if (results == NULL) {
    fprintf (stderr, "Tried to allocate %ld * %lld bytes for results_t and failed\n", sizeof(result_t), estimated_results);
    exit(222);
  }

  for (cur_phen=phenotypes; cur_phen != NULL; cur_phen = cur_phen->next) {
    fprintf (stderr, "Doing phenotype %s (iter %d)\n", cur_phen->name, iteration++);
    for (cur_snp = genotypes; cur_snp != NULL; cur_snp = cur_snp->next) {
      is_cis = check_cis(cur_snp, cur_phen, maxdist);
      if (is_cis == 0 && cis_only == 1) continue;
      switch (test_type) {
      case 0:
	p = nonparam_compar(cur_phen->values, cur_snp->gt, num_indivs, cur_snp->num_groups, sort_index, rank, tie_counts, &flag);
	break;
      case 1 :
	p = regression_significance (cur_snp->gt, cur_phen->values, num_indivs);
	flag = 0;
	break;
      default :
	Die("No such test type %d\n", test_type);
      }
      total_tests++;

      if (is_cis == 1) {
	total_cis_tests++;
      }
      if (p > 1.001) {
	fprintf (stderr, "P>1\trs%d\t%s\t%g\n", cur_snp->rs, cur_phen->name, p);
      }
      if (p <= MAXP) {
	if (tot_results < estimated_results) {
	  results[(tot_results)] = MallocOrDie(sizeof(result_t));
	  results[(tot_results)]->snp = cur_snp;
	  results[(tot_results)]->phen = cur_phen;
	  results[tot_results]->p = p;
	  results[tot_results]->flag = flag;
	  results[tot_results]->good_for_cis = is_cis;
	  tot_results++;
	} else {
	  Die("Estimated results off\n");
	}
      }
    }
  }

  *total_cis_tests_r = total_cis_tests;
  *tot_results_r = tot_results;
  *total_tests_r = total_tests;

  /* Now, sort the results in anticipation of B-H FDR */
  qsort (results, tot_results, sizeof(result_t *), &result_sort_func);

  if (test_type == 0) {
    free(sort_index);
    free(tie_counts);
    free(rank);
  }

  return(results);
}

void print_results (result_t **results, long long tot_results, double total_tests_d, double total_cis_tests_d, int cis_only) {

  long long fdr_threshold_index = -1;
  long long cis_fdr_threshold_index = -1;
  long long k_for_cis_fdr = 0;

  int i;
  int result_sig;

  /* Now, do B-H to find FDR threshold, both cis and trans */
  for (i=0; i < tot_results; i++) {
    if (results[i]->p <= ((double)(i+1.))/total_tests_d * FDR_ALPHA) {
      fdr_threshold_index = i;
    }
    if (results[i]->good_for_cis == 1) {
      k_for_cis_fdr++;
      if (results[i]->p <= ((double)(k_for_cis_fdr))/total_cis_tests_d * FDR_ALPHA) {
	cis_fdr_threshold_index = i;
      }
    }
  }

  /* Now, print results using bitwise marking of why sig:
     1 = Trans Bonferonni
     2 = Trans FDR
     4 = P<1e-05
     8 = Cis Bonferonni
     16 = Cis FDR */
  for (i = 0; i < tot_results; i++) {
    result_sig = 0;
    if (cis_only == 0) {
      if (results[i]->p < ALPHA/ total_tests_d) result_sig++;
      if (i <= fdr_threshold_index) result_sig += 2;
      if (results[i]->p < THRESHOLD) result_sig += 4;
    }
    if (results[i]->good_for_cis == 1) {
      if (results[i]->p < ALPHA/total_cis_tests_d) result_sig += 8;
      if (i <= cis_fdr_threshold_index) result_sig += 16;
    }
    if (result_sig > 0) {
      printf ("rs%d\t%d:%d\t%s\t%d:%d-%d\t%g\t%d\t%d\n", results[i]->snp->rs, results[i]->snp->chr, results[i]->snp->pos,
	      results[i]->phen->name, results[i]->phen->chr, results[i]->phen->start, results[i]->phen->stop,
	      results[i]->p, results[i]->flag, result_sig);
    }
  }
}

int main (int argc, char **argv) {
  snp_t *genotypes;
  phen_t *phenotypes;

  long long total_tests;
  long long total_cis_tests;

  result_t **results;
  long long tot_results;

  char *optname;                /* name of option found by Getopt()        */
  char *optarg;                 /* argument found by Getopt()              */
  int   optind;                 /* index in argv[]                         */
  
  int test_type = 0;            /* 0 = K-W, 1 = linear regression */
  int quant_norm = 0;           /* Do quantile normalization */
  int maxdist = 200000;         /* Max dist for cis search */
  int cis_only = 0;

  char *plink_prefix;
  char *gene_list;
  char *exp_dir;

  /**********************************************
   * Print header here 
   *********************************************/
  printf ("EQTL version %s\n", VERSION);
  printf ("%s\n%s\n\n", COPYRIGHT, LICENSE);

  /*********************************************** 
   * Parse command line
   ***********************************************/
  while (Getopt(argc, argv, OPTIONS, NOPTIONS, usage,
                &optind, &optname, &optarg))  {
    if (strcmp (optname, "--test") == 0) {
      if (strcmp(optarg, "kw") == 0) {
	test_type = 0;
      } else if (strcmp(optarg, "reg") == 0) {
	test_type = 1;
      } else {
	Die("Unrecognized test %s\n", optarg);
      }
    } else if (strcmp(optname, "--qnorm") == 0) {
      quant_norm = 1;
    } else if (strcmp(optname, "-c") == 0) {
      cis_only = 1;
    } else if (strcmp (optname, "--dist") == 0) {
      maxdist = 1000 * atoi(optarg);
    } else if (strcmp (optname, "-h") == 0) {
      puts(usage);
      puts(experts);
      exit(EXIT_SUCCESS);
    } 
  }
  if (argc - optind != 3) {
    Die("Incorrect number of arguments\n%s\n", usage);
  }

  plink_prefix=argv[optind++];
  gene_list =argv[optind++];
  exp_dir = argv[optind++];

  if (sizeof(long long) < 8) Die("Long is only %d; fix tot_tests\n", sizeof(long));

  genotypes = read_genotypes(plink_prefix);

  phenotypes = read_phenotypes (gene_list, exp_dir, genotypes->num_indivs, genotypes->id_list, quant_norm);

  results = get_results (genotypes, phenotypes, &tot_results, &total_tests, &total_cis_tests, test_type, cis_only, maxdist);

  printf ("There are %lld total tests and %lld total cis tests\n", total_tests, total_cis_tests);

  print_results (results, tot_results, (double)total_tests, (double)total_cis_tests, cis_only);
  
  printf ("\nFin\n");

  return(0);
}
  
