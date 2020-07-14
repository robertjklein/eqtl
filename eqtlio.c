/*
 * eqtlio.c
 *
 * Functions for reading in the files for eqtl analysis
 *
 * Robert J. Klein
 * July 8, 2010
 */

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <errno.h>
#include <gsl/gsl_cdf.h>

#include "squid.h"
#include "sqfuncs.h"

#include "structs.h"
#include "eqtlio.h"

char get_gt_code (char *c) {
  char *d;
  if (*c == '0') {
    return(0);
  } else {
    d = c;
    while (!isspace(*d)) d++;
    while (isspace(*d)) d++;
    if (*c == 'A')
      if (*d == 'A') return(1);
      else if (*d == 'C') return(2);
      else if (*d == 'G') return(3);
      else if (*d == 'T') return(4);
      else return (0);
    else if (*c == 'C') 
      if (*d == 'A') return(2);
      else if (*d == 'C') return(5);
      else if (*d == 'G') return(6);
      else if (*d == 'T') return(7);
      else return(0);
    else if (*c == 'G') 
      if (*d == 'A') return (3);
      else if (*d == 'C') return(6);
      else if (*d == 'G') return(8);
      else if (*d == 'T') return(9);
      else return(0);
    else if (*c == 'T')
      if (*d == 'A') return(4);
      else if (*d == 'C') return(7);
      else if (*d == 'G') return(9);
      else if (*d == 'T') return(10);
      else return(0);
    else return(0);
  }
}

void finish_recode (char *gt, int n, int g1, int g2, int g3) {
  int i, g;

  for (i=0; i<n; i++) {
    g = (int)(gt[i]);
    if (g == 0) {
      gt[i] = (char)127;
    } else if (g == g1) {
      gt[i] = (char)0;
    } else if (g == g2) {
      gt[i] = (char)1;
    } else if (g == g3) {
      gt[i] = (char)2;
    } else {
      Die("Bad char %c for %d,%d,%d\n", gt[i], g1,g2,g3);
    }
  }
}


/* Goal here is to take initial GT codings and recode for use in both K-W
   and regression.  Missing is coded as 127.
   Algorithm:
   1.  Count how many of each GT group.  
   2.  Create a bitmask to mark which groups are > 0
   3.  Manually code all valid 3- and 2-group mixtures.  Rest give errors
*/
int recode_gt (char *gt, int n) {
  int i;
  int counts[11];
  int cur_gt;
  int bitmask;

  for (i=0; i<12; i++) counts[i] = 0;
  for (i=0; i<n; i++) {
    cur_gt = (int)gt[i];
    if (cur_gt < 0 || cur_gt > 10) {
      Die("gt out of range (%d)\n", cur_gt);
    }
    counts[cur_gt]++;
  }

  /* Now, create bitmask */
  bitmask = 0;
  for (i=1; i<11; i++) {
    if (counts[i] > 0) {
      bitmask += (1 << (i-1));
    }
  }

  /* Now, manually test bitmask */
  switch (bitmask) {
    /* First, 3-group case */
  case (19):  /* A/C */
    finish_recode(gt, n, 1, 2, 5);
    return(3);
  case (133): /* A/G */
    finish_recode(gt, n, 1, 3, 8);
    return(3);
  case (521): /* A/T */
    finish_recode(gt, n, 1, 4, 10);
    return(3);
  case (176): /* C/G */
    finish_recode(gt, n, 5, 6, 8);
    return(3);
  case (592): /* C/T */
    finish_recode(gt, n, 5, 7, 10);
    return(3);
  case (896): /* G/T */
    finish_recode(gt, n, 8, 9, 10);
    return(3);
    
    /* Now, 2-group including first as het */
  case (3):   /* A/C */
    finish_recode (gt, n, 1, 2, 15);
    return(2);
  case (5): /* A/G */
    finish_recode (gt, n, 1, 3, 15);
    return(2);
  case (9): /* A/T */
    finish_recode (gt, n, 1, 4, 15);
    return(2);
  case (48): /* C/G */
    finish_recode (gt, n, 5, 6, 15);
    return(2);
  case (80): /* C/T */
    finish_recode (gt, n, 5, 7, 15);
    return(2);
  case (384): /* G/T */
    finish_recode (gt, n, 8, 9, 15);
    return(2);

    /* Now, 2-group including last as het */
  case (18):   /* A/C */
    finish_recode (gt, n, 2, 5, 15);
    return(2);
  case (132): /* A/G */
    finish_recode (gt, n, 3, 8, 15);
    return(2);
  case (520): /* A/T */
    finish_recode (gt, n, 4, 10, 15);
    return(2);
  case (160): /* C/G */
    finish_recode (gt, n, 6, 8, 15);
    return(2);
  case (576): /* C/T */
    finish_recode (gt, n, 7, 10, 15);
    return(2);
  case (768): /* G/T */
    finish_recode (gt, n, 9, 10, 15);
    return(2);
    
    /* Now, 2-groups, no hets */
  case (17):   /* A/C */
    finish_recode (gt, n, 1, 5, 15);
    return(2);
  case (129): /* A/G */
    finish_recode (gt, n, 1, 8, 15);
    return(2);
  case (513): /* A/T */
    finish_recode (gt, n, 1, 10, 15);
    return(2);
  case (144): /* C/G */
    finish_recode (gt, n, 5, 8, 15);
    return(2);
  case (528): /* C/T */
    finish_recode (gt, n, 5, 10, 15);
    return(2);
  case (640): /* G/T */
    finish_recode (gt, n, 8, 10, 15);
    return(2);
  }

  for (i=0; i<11; i++){
    printf ("%d\t%d\n", i, counts[i]);
  }
  Die("Bad bitmask %d\n", bitmask);
  return(0);
}

void check_recoded_gt (char *gt, int n, int g) {
  int i;
  int a;
  if (g == 2) {
    for (i=0; i<n; i++) {
      a = (int)gt[i];
      if (a != 0 && a != 1 && a != 15 && a != 127) {
	Die("Bad check\n");
      }
    }
  } else if (g == 3) {
    for (i=0; i<n; i++) {
      a = (int)gt[i];
      if (a != 0 && a != 1 && a != 2 && a != 127) {
	Die("Bad check\n");
      }
    }
  }
  printf ("%d\n", n);
}

snp_t *read_genotypes (char *filename) {
  snp_t *start, *cur, *prev;
  char **ids;
  int cur_id;
  FILE *f;
  char *buf;
  char *cp;
  int num_snps = 0, num_indivs = 0;

  buf = MallocOrDie(256);
  prev = NULL;

  sprintf (buf, "%s.map", filename);
  f = fopen(buf, "r");
  if (f==NULL) Die("Cannot open %s\n", buf);
  start = NULL;

  while (fgets(buf, 255, f)) {
    cur = (snp_t *)MallocOrDie(sizeof(snp_t));
    if (start == NULL) {
      start = cur;
    } else {
      prev->next = cur;
    }
    prev = cur;
    
    cur->next = NULL;
    cur->id_list = NULL;
    cur->gt = NULL;

    if ((isdigit(buf[0]) && isspace(buf[1])) ||
	(isdigit(buf[0]) && isdigit(buf[1]) && isspace(buf[2]))) {
      cur->chr = atoi(buf);
    } else if (buf[0] == 'X') {
      cur->chr = 23;
    } else if (buf[0] == 'Y') {
      cur->chr = 24;
    } else {
      cur->chr = 0;
    }
    
    cp = buf;
    while (!isspace(*cp)) cp++;
    while(isspace(*cp)) cp++;
    if (cp[0] == 'r' && cp[1] == 's') {
      cur->rs=atoi(cp+2);
    } else {
      cur->rs = 0;
    }
    
    while (!isspace(*cp)) cp++;
    while(isspace(*cp)) cp++;
    while (!isspace(*cp)) cp++;
    while (isspace(*cp)) cp++;
    cur->pos = atoi(cp);

    num_snps++;
  }
  fclose(f);

  free(buf);
  buf=MallocOrDie(sizeof(char)*(num_snps*4 + 256));

  sprintf (buf, "%s.ped", filename);
  f = fopen(buf, "r");
  if (f==NULL) Die("Cannot open %s\n", buf);

 

  while (fgets(buf, num_snps*4 + 254, f)) {
    for (cp=buf; *cp != 0 && *cp != '\n'; cp++);
    if (*cp == 0) {
      fprintf (stderr, "Length of %d\n", (int)(cp-buf));
      Die("Didn't read entire line\n");
    } else {
      /*      fprintf (stderr, "Did a line\n");*/
    }
    num_indivs++;
  }
  rewind(f);

  ids = MallocOrDie(sizeof(char *)*num_indivs);
  cur_id = 0;

  while (fgets(buf, num_snps*4 + 254, f)) {
    cp = buf;
    while (!isspace(*cp)) cp++;
    while (isspace(*cp)) cp++;
    while (!isspace(*cp)) cp++;
    ids[cur_id] = MallocOrDie((sizeof(char)*(cp-buf))+1);
    strncpy(ids[cur_id], buf, cp-buf);
    ids[cur_id][cp-buf] = '\0';

    while (isspace(*cp)) cp++;
    while (!isspace(*cp)) cp++;
    while (isspace(*cp)) cp++;
    while (!isspace(*cp)) cp++;
    while (isspace(*cp)) cp++;
    while (!isspace(*cp)) cp++;
    while (isspace(*cp)) cp++;
    while (!isspace(*cp)) cp++;
    while (isspace(*cp)) cp++;

    /* Now, we're at the gt's */
    for (cur=start; cur != NULL; cur = cur->next) {
      if (cur->gt == NULL) {
	cur->gt = MallocOrDie(sizeof(char)*num_indivs);
      }
      cur->gt[cur_id] = get_gt_code(cp);
      while (!isspace(*cp)) cp++;
      while (isspace(*cp)) cp++;
      while (!isspace(*cp)) cp++;
      while (isspace(*cp)) cp++;
    }

    cur_id++;
  }


  for (cur=start; cur!=NULL; cur = cur->next) {
    cur->num_indivs = num_indivs;
    cur->num_snps = num_snps;
    cur->id_list = ids;
    /*printf ("Recoding a gt rs%d\n", cur->rs);*/
    cur->num_groups = recode_gt (cur->gt, num_indivs);
    /*check_recoded_gt (cur->gt, num_indivs, cur->num_groups);*/
  }
  free(buf);
  return(start);
}

int val_sort_func (const void *a, const void *b) {
  float i,j;

  i = **(float **)a;
  j = **(float **)b;
  
  if (i < j) {
    return(-1);
  } else if (i > j) {
    return(1);
  } else {
    return(0);
  }
}

void quantile_normalize (float *vals, int n) {
  float **sort_index;
  int i;
  float increment;
  
  sort_index = MallocOrDie(sizeof(float *)*n);
  for (i=0; i<n; i++) sort_index[i] = vals+i;

  qsort (sort_index, n, sizeof(float *), &val_sort_func);

  increment = 1./n;

  for (i=0; i<n; i++) {
    *(sort_index[i]) = (float)gsl_cdf_ugaussian_Pinv((double)((i+0.5)*increment));
  }

}

phen_t *read_phenotypes (char *probelist, char *probedir, int num_indivs, char **id_list, int qnorm) {
  char buf[256];
  FILE *f;
  phen_t *start, *cur, *prev;
  char *cp;
  int i;
  int tot_read;

  start = NULL;
  prev = NULL;

  f = fopen(probelist, "r");
  if (f==NULL) Die("Cannot open %s\n", probelist);
  
  while (fgets(buf, 255, f)) {
    cur = MallocOrDie(sizeof(phen_t));
    if (start == NULL) {
      start = cur;
    } else {
      prev->next = cur;
    }
    prev = cur;

    /* Initial copy of name */
    for (cp=buf; !isspace(*cp); cp++);
    *cp = '\0';
    cur->name = MallocOrDie(sizeof(char)*(strlen(buf)+1));
    strncpy(cur->name, buf, strlen(buf)+1);
 
    /* Set chr, start stop */
    cp++;
    while (isspace(*cp)) cp++;
    if (*cp == 'X') 
      cur->chr = 23;
    else
      cur->chr = atoi(cp);
    while (!isspace(*cp)) cp++;
    while (isspace(*cp)) cp++;
    cur->start = atoi(cp);
    while (!isspace(*cp)) cp++;
    while (isspace(*cp)) cp++;
    cur->stop = atoi(cp);

    for (cp=cur->name; *cp != '\0' && !isspace(*cp) && isprint (*cp); cp++);
    *cp = '\0';
    
    cur->values = MallocOrDie (sizeof(float)*num_indivs);
    for (i=0; i<num_indivs; i++) {
      cur->values[i] = 0;
    }
    cur->next = NULL;
  }
  fclose(f);

  for (cur=start; cur != NULL; cur = cur->next) {
    sprintf (buf, "%s/%s.phen", probedir, cur->name);
    f = fopen(buf, "r");
    if (f==NULL) {
      for (cp=buf; *cp != '\0'; cp++) {
	fprintf (stderr, "%c %d\n", *cp, (int)(*cp));
      }
      Die("Could not open %s b/c of %d\n", buf, errno);
    }
    tot_read = 0;
    while (fgets (buf, 255, f)) {
      for (i=0; i < num_indivs && strncmp(id_list[i], buf, strlen(id_list[i])) != 0; i++);
      if (i < num_indivs) {
	cp = buf + strlen(id_list[i]);
	while (isspace(*cp)) cp++;
	if (isdigit(*cp) || *cp == '-') {
	  cur->values[i] = atof(cp);
	} else {
	  fprintf (stderr, "WARNING: Found a non-number file %s/%s.phen line %s", probedir, cur->name, buf);
	  cur->values[i] = 0.;
	}
	tot_read++;
      } 
    }
    fclose(f);
    if (tot_read < num_indivs) Die("Not enough individuals in %s\n", cur->name);
    if (qnorm == 1) quantile_normalize(cur->values, num_indivs);
  }
  return(start);
}
