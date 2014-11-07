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

#ifndef SRC_TESTS_LOGISTIC_REGRESSION_TEST_H__
#define SRC_TESTS_LOGISTIC_REGRESSION_TEST_H__

#include <cppunit/extensions/HelperMacros.h>

#include "src/logistic_regression.h"

class LogisticRegressionTest :
public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(LogisticRegressionTest);
  CPPUNIT_TEST(test_logistic);
  CPPUNIT_TEST(test_logistic_regression);
  CPPUNIT_TEST(test_logistic_regression_prediction);
  CPPUNIT_TEST_SUITE_END();

 public:
  void setUp();
  void tearDown();
 private:
  void test_logistic();
  void test_logistic_regression();
  void test_logistic_regression_prediction();
};

#endif //  SRC_TESTS_LOGISTIC_REGRESSION_TEST_H__
