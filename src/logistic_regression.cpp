/*
Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>

This file is part of lobSTR.

lobSTR is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

lobSTR is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with lobSTR.  If not, see <http://www.gnu.org/licenses/>.

*/

#include <math.h>
#include <gsl/gsl_multimin.h>

#include "src/logistic_regression.h"

using namespace std;

double logistic(const double z) {
  return 1.0/(1+exp(-1*z));
}

double logistic_regression_f(const gsl_vector *theta, void *params) {
  // Pull out data from params
  double *data = (double *)params;
  int m = data[0];
  int n = data[1];
  // Get y
  double y[m];
  for (int i = 0; i < m; i++) {
    y[i] = data[i+2];
  }
  // Get theta
  double coeffs[n+1];
  for (int i = 0; i < n+1; i++) {
    coeffs[i] = gsl_vector_get(theta, i);
  }
  // compute cost
  double cost = 0;
  for (int i = 0; i < m; i++) { // i counts which data point we're on
    double z = coeffs[0];
    for (int j = 0; j < n; j++) { // j counts which param
      z += data[m+i*n+j+2]*coeffs[j+1];
    }
    double h = logistic(z);
    cost += 1/(float)m*(-1*y[i]*log(h)-(1-y[i])*log(1-h));
  }
  return cost;
}

void logistic_regression_df(const gsl_vector *theta, void *params, gsl_vector *df) {
  // Pull out data from params
  double *data = (double *)params;
  int m = data[0];
  int n = data[1];
  // Get y
  double y[m];
  for (int i = 0; i < m; i++) {
    y[i] = data[i+2];
  }
  // Get theta
  double coeffs[n+1];
  for (int i = 0; i < n+1; i++) {
    coeffs[i] = gsl_vector_get(theta, i);
  }
  // Init gradient to 0
  for (int k = 0; k < n+1; k++) {
    gsl_vector_set(df, k, 0);
  }
  // Update
  for (int i = 0; i < m; i++) { // i counts which data point we're on
    double z = coeffs[0];
    for (int j = 0; j < n; j++) { // j counts which param
      z += data[m+i*n+j+2]*coeffs[j+1];      
    }
    double h = logistic(z);
    double err = h-y[i];
    // Update each dimension for this data point
    for (int k = 0; k < n+1; k++) {
      double xval = 1;
      if (k != 0) {
	xval = data[2+m+i*n+(k-1)];
      }
      gsl_vector_set(df, k, gsl_vector_get(df, k)+err*xval*1/(float)m);
    }
  }
}

void logistic_regression_fdf(const gsl_vector *x, void* params,
	    double *f, gsl_vector *df) {
  *f = logistic_regression_f(x, params);
  logistic_regression_df(x, params, df);
}

int logistic_regression(const vector<bool>& y, const vector<vector<double> >& x,
			const int& m, const int& n,
			const double& step_size, const double& tol,
			const double& epsabs, const size_t& maxiter,
			vector<double>* coeffs) {
  int retcode = 0;
  // Set up parameters
  gsl_vector *theta = gsl_vector_alloc(n+1);
  gsl_vector_set_zero(theta);
  // Set up data
  vector<double> data (2+m*(n+1), 0);
  data[0] = m;
  data[1] = n;
  copy(y.begin(), y.end(), data.begin()+2);
  for(size_t i=0; i<x.size(); i++) {
    copy(x[i].begin(), x[i].end(), data.begin()+2+m+i*n);
  }
  // Set up minimizer
  size_t iter = 0;
  int status;
  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;
  gsl_multimin_function_fdf logreg;
  logreg.f = &logistic_regression_f;
  logreg.df = &logistic_regression_df;
  logreg.fdf = &logistic_regression_fdf;
  logreg.n = n+1;
  logreg.params = &data[0];
  T = gsl_multimin_fdfminimizer_conjugate_fr;
  s = gsl_multimin_fdfminimizer_alloc(T, n+1);
  gsl_multimin_fdfminimizer_set(s, &logreg, theta, step_size, tol);
  // Run the minimizer
  do {
    iter++;
    status = gsl_multimin_fdfminimizer_iterate(s);
    if (status) {
      retcode = 1;
      break;
    }
    status = gsl_multimin_test_gradient(s->gradient, epsabs);
  } while (status == GSL_CONTINUE && iter < maxiter);
  if (iter >= maxiter) retcode = 1;
  for (int k = 0; k < n+1; k++) {
    coeffs->at((size_t)k) = gsl_vector_get(s->x, k);
  }
  gsl_vector_free(theta);
  return retcode;
}

double logistic_regression_predict(std::vector<double> x, std::vector<double> coeffs) {
  size_t n = coeffs.size();
  double z = coeffs.at(0);
  for (size_t i=0; i<n-1; i++) {
    z += x.at(i)*coeffs.at(i+1);
  }
  return logistic(z);
}
