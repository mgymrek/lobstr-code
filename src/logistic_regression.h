/*
Copyright (C) 2011-2014 Melissa Gymrek <mgymrek@mit.edu>

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

#include <gsl/gsl_vector.h>

#include <vector>

#ifndef SRC_LOGISTIC_REGRESSION_H_
#define SRC_LOGISTIC_REGRESSION_H_

/*
Compute the logistic function
 */
double logistic(const double z);

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
double logistic_regression_f(const gsl_vector *theta, void *params);

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
void logistic_regression_df(const gsl_vector *theta, void *params, gsl_vector *df);

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
			     double *f, gsl_vector *df);


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
int logistic_regression(const std::vector<bool>& y, const std::vector<std::vector<double> >& x,
			const int& m, const int& n,
			const double& step_size, const double& tol,
			const double& epsabs, const size_t& maxiter,
			std::vector<double>* coeffs);

/*
Predict probability using logistic regression model

Inputs:
vector<double> x: independent variables
vector<double> coeffs: coefficients
 */
double logistic_regression_predict(std::vector<double> x, std::vector<double> coeffs);

#endif  // SRC_LOGISTIC_REGRESSION_H_
