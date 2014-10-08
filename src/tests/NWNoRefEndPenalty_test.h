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

#ifndef SRC_TESTS_NWNOREFENDPENALTY_H__
#define SRC_TESTS_NWNOREFENDPENALTY_H__

#include <cppunit/extensions/HelperMacros.h>

#include "src/NWNoRefEndPenalty.h"

class NWNoRefEndPenaltyTest :
public CppUnit::TestFixture {
  CPPUNIT_TEST_SUITE(NWNoRefEndPenaltyTest);
  CPPUNIT_TEST(test_Align_Mut_Only);
  CPPUNIT_TEST(test_Align_Ins_Only);
  CPPUNIT_TEST(test_Align_Del_Only);
  CPPUNIT_TEST(test_Align_All_Muts);
  CPPUNIT_TEST_SUITE_END();

 public:
  static const int MAX_INS=10, MAX_DEL=10, PERFECT_FLANK=10; 
  static const size_t REF_LEN=150, READ_LEN=50; 
  static const int RAND_SEED=21312278;
  static const int NUM_TRIALS=1000;
  static const double MIN_FRAC_CORRECT=0.9;

  void setUp();
  void tearDown();
  int  GenAlignments(int num_trials, double mut_prob, double ins_prob, double del_prob);
  void test_Align_Mut_Only();
  void test_Align_Ins_Only();
  void test_Align_Del_Only();
  void test_Align_All_Muts();
  
 private:
};

#endif //  SRC_TESTS_NWNOREFENDPENALTY_H_
