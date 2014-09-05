#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <math.h>
#include <err.h>
#include <iostream>
#include <vector>
#include "src/logistic_regression.h"

using namespace std;

/*
Should match this code in R

Example 1 (1 ind. var)
x = c(-1,5,3,-2,3,4,1,2,3,1)
y = c(0,0,0,1,0,1,0,1,0,1)
log.model = glm(y~x, family="binomial")
summary(log.model)$coefficients[c(1,2)]
[1]  0.06485105 -0.25951130

Example 2 (2 ind. var)
x1 = c(-1,5,3,-2,3,4,1,2,3,1)
x2 = c(-10,4,2,-4,3,5,6,1,3,1)
y = c(0,0,0,1,0,1,0,1,0,1)
log.model = glm(y~x1+x2, family="binomial")
summary(log.model)$coefficients[c(1,2,3)]
[1]  0.3459375 -0.4916129  0.1354238
 */
int main (void)
{
  // Logistic regression params
  double step_size = 0.001;
  double tol = 1e-4;
  double epsabs = 1e-3;
  size_t maxiter = 100;

  // Example 1 - 1 independent variable
  int m = 10;
  int n = 1;
  bool yy[] = {0,0,0,1,0,1,0,1,0,1};
  double xx[] = {-1,5,3,-2,3,4,1,2,3,1};
  vector<bool> y;
  y.assign(yy, yy+m);
  vector<vector <double> >x;
  for (int i=0; i<m; i++) {
    vector<double> x0;
    x0.push_back(xx[i]);
    x.push_back(x0);
  }
  vector<double> coeffs;
  coeffs.resize(n+1);
  logistic_regression(y, x, m, n, step_size, tol, epsabs, maxiter, &coeffs);
  printf("Coefficient results: b0=%g, b1=%g (R results: 0.06485105, -0.25951130\n", coeffs.at(0), coeffs.at(1));
  if (coeffs.at(0) < 0.0647 || coeffs.at(1) > 0.0649) {
    errx(1, "Example 1: Intercept term doesn't match");
  }
  if (coeffs.at(1) < -0.26 || coeffs.at(1) > -0.258) {
    errx(1, "Example 1: Slope term doesn't match");
  }

  // Example 2 - 2 independent variables
  m = 10;
  n = 2;
  double xx1[] = {-1,5,3,-2,3,4,1,2,3,1};
  double xx2[] = {-10,4,2,-4,3,5,6,1,3,1};
  x.clear();
  for (int i=0; i<m; i++) {
    vector<double>x0;
    x0.push_back(xx1[i]);
    x0.push_back(xx2[i]);
    x.push_back(x0);
  }
  coeffs.resize(n+1);
  logistic_regression(y, x, m, n, step_size, tol, epsabs, maxiter, &coeffs);
  printf("Coefficient results: b0=%g, b1=%g, b2=%g (R results: 0.3459375, -0.4916129, 0.1354238\n", coeffs.at(0), coeffs.at(1), coeffs.at(2));
  if (coeffs.at(0) < 0.33 || coeffs.at(0) > 0.35) {
    errx(1, "Example 2: Intercept term doesn't match");
  }
  if (coeffs.at(1) < -0.50 || coeffs.at(1) > -0.48) {
    errx(1, "Example 2: Slope1 term doesn't match");
  }
  if (coeffs.at(2) < 0.12 || coeffs.at(2) > 0.14) {
    errx(1, "Example 2: Slope2 term doesn't match");
  }
  return 0;
}
