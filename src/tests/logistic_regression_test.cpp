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

#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "src/tests/logistic_regression_test.h"
#include "src/runtime_parameters.h"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(LogisticRegressionTest);

// Logistic regression params
const double step_size = 0.001;
const double tol = 1e-4;
const double epsabs = 1e-3;
const size_t maxiter = 100;

void LogisticRegressionTest::setUp() {}

void LogisticRegressionTest::tearDown() {}

void LogisticRegressionTest::test_logistic() {
  CPPUNIT_ASSERT_MESSAGE("logistic failed", abs(logistic(2.3)-0.908877) < 0.01);
  CPPUNIT_ASSERT_MESSAGE("logistic failed", abs(logistic(0)-0.5) < 0.01);
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
void LogisticRegressionTest::test_logistic_regression() {
  // Example 1 - 1 independent variable
  int m = 10;
  int n = 1;
  bool yy[] = {0,0,0,1,0,1,0,1,0,1};
  double xx[] = {-1,5,3,-2,3,4,1,2,3,1};
  std::vector<bool> y;
  y.assign(yy, yy+m);
  std::vector<std::vector <double> >x;
  for (int i=0; i<m; i++) {
    std::vector<double> x0;
    x0.push_back(xx[i]);
    x.push_back(x0);
  }
  std::vector<double> coeffs;
  coeffs.resize(n+1);
  logistic_regression(y, x, m, n, step_size, tol, epsabs, maxiter, &coeffs);
  CPPUNIT_ASSERT_MESSAGE("Example 1: Intercept term doesn't match", abs(coeffs.at(0)-0.06485)<0.01);
  CPPUNIT_ASSERT_MESSAGE("Example 1: Slope term doesn't match", abs(coeffs.at(1)-(-0.257))<0.01);

  // Example 2 - 2 independent variables
  m = 10;
  n = 2;
  double xx1[] = {-1,5,3,-2,3,4,1,2,3,1};
  double xx2[] = {-10,4,2,-4,3,5,6,1,3,1};
  x.clear();
  for (int i=0; i<m; i++) {
    std::vector<double>x0;
    x0.push_back(xx1[i]);
    x0.push_back(xx2[i]);
    x.push_back(x0);
  }
  coeffs.resize(n+1);
  logistic_regression(y, x, m, n, step_size, tol, epsabs, maxiter, &coeffs);
  CPPUNIT_ASSERT_MESSAGE("Example 2: Intercept term doesn't match", abs(coeffs.at(0)-0.34)<0.01);
  CPPUNIT_ASSERT_MESSAGE("Example 2: Slope1 term doesn't match", abs(coeffs.at(1)-0.49)<0.01);
  CPPUNIT_ASSERT_MESSAGE("Example 2: Slope2 term doesn't match", abs(coeffs.at(2)-0.13)<0.01);
}

void LogisticRegressionTest::test_logistic_regression_prediction() {
  std::vector<double> vals;
  vals.push_back(1);
  vals.push_back(2);
  std::vector<double> cfs;
  cfs.push_back(-0.5);
  cfs.push_back(1.2);
  cfs.push_back(3);
  double pred = logistic_regression_predict(vals, cfs);
  CPPUNIT_ASSERT_MESSAGE("Logistic regression prediction failed", abs(pred-0.99877)<0.01);
}



