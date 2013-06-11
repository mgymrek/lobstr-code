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

#ifndef SRC_TESTS_COMMON_TEST_H__
#define SRC_TESTS_COMMON_TEST_H__

#include <cppunit/extensions/HelperMacros.h>

#include "src/common.h"

class CommonTest :
public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(CommonTest);
  CPPUNIT_TEST(test_TrimRead);
  CPPUNIT_TEST(test_getMSSeq);
  CPPUNIT_TEST(test_getCanonicalRepeat);
  CPPUNIT_TEST(test_getCanonicalMS);
  CPPUNIT_TEST(test_IsPerfectRepeat);
  CPPUNIT_TEST(test_reverseComplement);
  CPPUNIT_TEST(test_reverse);
  CPPUNIT_TEST_SUITE_END();

 public:
  void setUp();
  void tearDown();
 private:
  void test_TrimRead();
  void test_getMSSeq();
  void test_getCanonicalRepeat();
  void test_getCanonicalMS();
  void test_IsPerfectRepeat();
  void test_reverseComplement();
  void test_reverse();
};

#endif //  SRC_TESTS_COMMON_TEST_H__
