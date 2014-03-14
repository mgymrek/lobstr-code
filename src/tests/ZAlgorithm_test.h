/*
Copyright (C) 2014 Thomas Willems <twillems@mit.edu>

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

#ifndef SRC_TESTS_ZALGORITHM_H__
#define SRC_TESTS_ZALGORITHM_H__

#include <cppunit/extensions/HelperMacros.h>

#include "src/ZAlgorithm.h"

class ZAlgorithmTest :
public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(ZAlgorithmTest);
  CPPUNIT_TEST(test_Prefix);
  CPPUNIT_TEST(test_Suffix);
  CPPUNIT_TEST_SUITE_END();

 public:
  const static int NUM_TRIALS=1000;
  const static int LENGTH=100;
  void setUp();
  void tearDown();
  void test_Prefix();
  void test_Suffix();
    
 private:
};

#endif //  SRC_TESTS_ZALGORITHM_H_
