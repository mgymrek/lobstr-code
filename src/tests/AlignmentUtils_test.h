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

#ifndef SRC_TESTS_ALIGNMENTUTILS_H__
#define SRC_TESTS_ALIGNMENTUTILS_H__

#include <cppunit/extensions/HelperMacros.h>

#include "src/AlignmentUtils.h"

class AlignmentUtilsTest :
public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(AlignmentUtilsTest);
  CPPUNIT_TEST(test_StitchReads);
  CPPUNIT_TEST(test_GetMapq);
  CPPUNIT_TEST(test_GetSTRAllele);
  CPPUNIT_TEST_SUITE_END();

 public:
  void setUp();
  void tearDown();

  void test_StitchReads();
  void test_GetMapq();
  void test_GetSTRAllele();
};

#endif //  SRC_TESTS_ALIGNMENTUTILS_H__
