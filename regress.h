/*
 * regress.h
 *
 * Interface to simple linear regression significance testing
 *
 * Robert J. Klein
 * July 25, 2007
 */

#ifndef _regress_h
#define _regress_h

float regression_significance (char *gts, float *vals, int n);

#endif
