/* 
 * regress.c 
 * 
 * Simple linear regression significance testing
 * 
 * Robert J. Klein
 * July 25, 2007
 *
 */

#include <stdlib.h>
#include <math.h>

#include <gsl/gsl_cdf.h>

/*
 * Computes regression line, and t-test for slope of line != 0
 * Algorithm taken from sections 7.8 and 8.8 of Hogg and Tanis,
 * Probability and Statistical Inference, 6th edition, and
 * validated by comparing on their test data with R's results on same
 * data
 */
float regression_significance (char *gts, float *vals, int n_tot) {
  double sum_xy, sum_y, sum_y2;
  int sum_x, sum_x2;
  double sum_x_xbar;
  int i;
  double beta_hat, n_sigma2_hat2;
  double mean;
  double t1;
  int n, g;
  double v;

  sum_xy = 0.;
  sum_y = 0.;
  sum_y2 = 0.;
  sum_x = 0;
  sum_x2 = 0;

  n = 0;

  for (i=0; i<n_tot; i++) {
    g = (int)(gts[i]);
    if (g != 127) {
      n++;

      sum_x += g;
      sum_x2 += g*g;

      v = (double)vals[i];
      sum_y += v;
      sum_y2 += v*v;
      sum_xy += v*g;
    }
  }
  beta_hat = (sum_xy-sum_x*(sum_y/n))/(sum_x2-(1./n)*sum_x*sum_x);

  n_sigma2_hat2 = sum_y2 - sum_y*sum_y/n - beta_hat*sum_xy + (beta_hat*sum_x)*(sum_y/n);

  sum_x_xbar = 0.;
  mean = (1./n)*sum_x;
  for (i=0; i<n_tot; i++) {
    g = (int)(gts[i]);
    if (g != 127) {
      sum_x_xbar += (g-mean)*(g-mean);
    }
  }

  t1 = beta_hat/sqrt(n_sigma2_hat2/((n-2)*sum_x_xbar));
  return((float)(2*gsl_cdf_tdist_Q((double)fabs(t1), (double)(n-2))));
  //return(1.);
}
  
  
  
