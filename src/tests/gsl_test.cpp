#include <stdio.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multimin.h>
#include <math.h>
#include <err.h>
#include <iostream>
#include <vector>

using namespace std;

/*
Compute the logistic function
 */
double logistic(const double z) {
  return 1.0/(1+exp(-1*z));
}

/*
Logistic regression cost function

J(theta) = mean(-y_i*log(h_i) - (1-y_i)*log(1-h_i))
where: h_i = logistic(theta*x_i)

Inputs:
gsl_vector *theta: consists of param values theta. length n+1
void *params: length 2+m*(n+1). consists of the following unrolled to a single vector:
  - m: the number of data points
  - n: the number of independent variables
  - y1...ym: the classification values. Must be 0 or 1
  - x_(1,1), x(1,2),...x(1,n), x(2,1), x(2,2)...x(2,n),...x(m,n): all m x-vectors, each of length n, concatenated

Outputs:
double cost: evaluation of the cost function at the given parameters
 */
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

/*
Logistic regression gradient function

d(J(theta))/d(theta_j) = 1/m*sum_i=1^m (h_i-y_i)*x_(i,j)
where: h_i = logistic(theta*x_i)

Inputs:
gsl_vector *theta: consists of param values theta. length n+1
void *params: length 2+m*(n+1). consists of the following unrolled to a single vector:
  - m: the number of data points
  - n: the number of independent variables
  - y1...ym: the classification values. Must be 0 or 1
  - x_(1,1), x(1,2),...x(1,n), x(2,1), x(2,2)...x(2,n),...x(m,n): all m x-vectors, each of length n, concatenated

Outputs:
gsl_vector df: evaluation of gradient wrt each coefficient. length n+1
 */
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

/*
  Compute both f and df together

Inputs:
gsl_vector *x: consists of param values theta. length n+1
void *params: length 2+m*(n+1). consists of the following unrolled to a single vector:
  - m: the number of data points
  - n: the number of independent variables
  - y1...ym: the classification values. Must be 0 or 1
  - x_(1,1), x(1,2),...x(1,n), x(2,1), x(2,2)...x(2,n),...x(m,n): all m x-vectors, each of length n, concatenated

Outputs:
double *f: evaluation of cost function
gsl_vector *df: evaluation of gradient
 */
void logistic_regression_fdf(const gsl_vector *x, void* params,
	    double *f, gsl_vector *df) {
  *f = logistic_regression_f(x, params);
  logistic_regression_df(x, params, df);
}

/*
Run logistic regression

Inputs;
vector<bool> y: binary vector of response variables
vector<vector<float>>x: atrix of independent variables
int m: number of data points
int n: number of independent variables
double step_size: first step size for following the gradient
double tol: accuracy of the line minimization (make search direction orthogonal to gradient to some tolerance tol)
double epsabs: return success if the gradient is smaller than this number

Outputs:
vector<double>* coeffs: vector of ML coefficients

Returns 0 if completed successfully, else 1.
 */
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
