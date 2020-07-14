/*
 * eqtlio.h
 *
 * IO funcs
 */

#ifndef _eqtlio_h
#define _eqtlio_h

#include "structs.h"

snp_t *read_genotypes (char *filename);

phen_t *read_phenotypes (char *probelist, char *probedir, int num_indivs, char **id_list, int qnorm);

#endif
