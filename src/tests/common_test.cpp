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

#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
#include <map>
#include <string>
#include <vector>

#include "src/tests/common_test.h"
#include "src/runtime_parameters.h"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(CommonTest);

void CommonTest::setUp() {}

void CommonTest::tearDown() {}

void CommonTest::test_TrimRead() {
  std::string input_nucs =  "ACAGTCGATCGTAGCTAGCGCTAGCTACTAGCATCGATCGATCGATCGTACGTACGTAGGCATCATCG";
  std::string input_quals = "eeeeeeeeeeeeeffffffffffffffffffffffffffffffffffffffff'''''''''''''''";
  std::string trimmed_nucs;
  std::string trimmed_quals;
  TrimRead(input_nucs, input_quals, &trimmed_nucs, &trimmed_quals, QUAL_CUTOFF);
  CPPUNIT_ASSERT_MESSAGE("wrong trimmed nucs", trimmed_nucs == "ACAGTCGATCGTAGCTAGCGCTAGCTACTAGCATCGATCGATCGATCGTACGT");
  CPPUNIT_ASSERT_MESSAGE("wrong trimmed quals", trimmed_quals == "eeeeeeeeeeeeeffffffffffffffffffffffffffffffffffffffff");
}

void CommonTest::test_CheckRepeatCount() {
  std::string bestkmer;
  bool check;
  check = CheckRepeatCount("ACTAGCTACTACGTACGTAGCTGA", 1, 10, &bestkmer);
  CPPUNIT_ASSERT_MESSAGE("CheckRepeatCount failed", bestkmer == "A");
  CPPUNIT_ASSERT_MESSAGE("CheckRepeatCount failed", !check);
  check = CheckRepeatCount("ACACACACACACACAC", 2, 10, &bestkmer);
  CPPUNIT_ASSERT_MESSAGE("CheckRepeatCount failed", bestkmer == "AC");
  CPPUNIT_ASSERT_MESSAGE("CheckRepeatCount failed", check);
  check = CheckRepeatCount("ACACACACACACACAC", 7, 10, &bestkmer);
  CPPUNIT_ASSERT_MESSAGE("CheckRepeatCount failed", !check);
  check = CheckRepeatCount("ACACACACACACACACA", 1, 10, &bestkmer);
  CPPUNIT_ASSERT_MESSAGE("CheckRepeatCount failed", bestkmer == "A");
  CPPUNIT_ASSERT_MESSAGE("CheckRepeatCount failed", !check);
}

void CommonTest::test_reverseComplement() {
  // Case 1: upper case
  std::string nucs = "ACGATCGTGTCATGCNNACCACG";
  std::string rev = reverseComplement(nucs);
  CPPUNIT_ASSERT_MESSAGE("Incorrect reverse complement", rev == "CGTGGTNNGCATGACACGATCGT");

  // Case 2: lower case
  nucs = "acgaccacagctacgacnacgactan";
  rev = reverseComplement(nucs);
  CPPUNIT_ASSERT_MESSAGE("Incorrect reverse complement", rev == "NTAGTCGTNGTCGTAGCTGTGGTCGT");
}

void CommonTest::test_reverse() {
  std::string quals = "ffbcefdcdffcfdffddaddadaeffc";
  std::string rev = reverse(quals);
  CPPUNIT_ASSERT_MESSAGE("Incorrect reverse", rev == "cffeadaddaddffdfcffdcdfecbff");
}

void CommonTest::test_ExtractCigar() {
  // Simple cigar: 50M5D50M. Test STR overlapping the deletion, on one end, or ending at boundary
  CIGAR_LIST cigar_list, str_cigar_list;
  cigar_list.cigars.clear();
  CIGAR cig;
  cig.num = 50;
  cig.cigar_type = 'M';
  cigar_list.cigars.push_back(cig);
  cig.num = 5;
  cig.cigar_type = 'D';
  cigar_list.cigars.push_back(cig);
  cig.num = 50;
  cig.cigar_type = 'M';
  cigar_list.cigars.push_back(cig);
  cigar_list.ResetString();
  ExtractCigar(cigar_list, 0, 40, 60, &str_cigar_list); // overlapping deletion
  std::string expected_cigar_string = "50M5D50M";
  CPPUNIT_ASSERT_MESSAGE("ExtractCigar failed", str_cigar_list.cigar_string == expected_cigar_string);
  ExtractCigar(cigar_list, 0, 10, 20, &str_cigar_list);
  expected_cigar_string = "50M"; // left flank
  CPPUNIT_ASSERT_MESSAGE("ExtractCigar failed", str_cigar_list.cigar_string == expected_cigar_string);
  ExtractCigar(cigar_list, 0, 70, 80, &str_cigar_list);
  expected_cigar_string = "50M"; // right flank
  CPPUNIT_ASSERT_MESSAGE("ExtractCigar failed", str_cigar_list.cigar_string == expected_cigar_string);
  ExtractCigar(cigar_list, 0, 50, 60, &str_cigar_list); // boundary
  expected_cigar_string = "50M5D50M";
  CPPUNIT_ASSERT_MESSAGE("ExtractCigar failed", str_cigar_list.cigar_string == expected_cigar_string);
  ExtractCigar(cigar_list, 100, 150, 160, &str_cigar_list); // non-zero start coord
  expected_cigar_string = "50M5D50M";
  CPPUNIT_ASSERT_MESSAGE("ExtractCigar failed", str_cigar_list.cigar_string == expected_cigar_string);
  
  // Two different indels: 50M2I60M5I90M
  cigar_list.cigars.clear();
  cig.num = 50;
  cig.cigar_type = 'M';
  cigar_list.cigars.push_back(cig);
  cig.num = 2;
  cig.cigar_type = 'I';
  cigar_list.cigars.push_back(cig);
  cig.num = 60;
  cig.cigar_type = 'M';
  cigar_list.cigars.push_back(cig);
  cig.num = 5;
  cig.cigar_type = 'I';
  cigar_list.cigars.push_back(cig);
  cig.num = 90;
  cig.cigar_type = 'M';
  cigar_list.cigars.push_back(cig);
  cigar_list.ResetString();
  ExtractCigar(cigar_list, 0, 30, 60, &str_cigar_list); // overlapping first insertion
  expected_cigar_string = "50M2I60M";
  CPPUNIT_ASSERT_MESSAGE("ExtractCigar failed", str_cigar_list.cigar_string == expected_cigar_string);
  ExtractCigar(cigar_list, 0, 100, 130, &str_cigar_list); // overlapping second insertion
  expected_cigar_string = "60M5I90M";
  CPPUNIT_ASSERT_MESSAGE("ExtractCigar failed", str_cigar_list.cigar_string == expected_cigar_string);
  ExtractCigar(cigar_list, 0, 150, 190, &str_cigar_list); // all in right flank
  expected_cigar_string = "90M";
  CPPUNIT_ASSERT_MESSAGE("ExtractCigar failed", str_cigar_list.cigar_string == expected_cigar_string);
  ExtractCigar(cigar_list, 0, 30, 190, &str_cigar_list); // overlapping both
  expected_cigar_string = "50M2I60M5I90M";
  CPPUNIT_ASSERT_MESSAGE("ExtractCigar failed", str_cigar_list.cigar_string == expected_cigar_string);
  ExtractCigar(cigar_list, 0, 110, 190, &str_cigar_list); // boundary
  expected_cigar_string = "60M5I90M";
  CPPUNIT_ASSERT_MESSAGE("ExtractCigar failed", str_cigar_list.cigar_string == expected_cigar_string);
  ExtractCigar(cigar_list, 0, 30, 50, &str_cigar_list); // other boundary
  expected_cigar_string = "50M2I60M";
  CPPUNIT_ASSERT_MESSAGE("ExtractCigar failed", str_cigar_list.cigar_string == expected_cigar_string);

  // Two different indels with indels at the ends 2I50M2I60M5I90M2I
  cigar_list.cigars.clear();
  cig.num = 2;
  cig.cigar_type = 'I';
  cigar_list.cigars.push_back(cig);
  cig.num = 50;
  cig.cigar_type = 'M';
  cigar_list.cigars.push_back(cig);
  cig.num = 2;
  cig.cigar_type = 'I';
  cigar_list.cigars.push_back(cig);
  cig.num = 60;
  cig.cigar_type = 'M';
  cigar_list.cigars.push_back(cig);
  cig.num = 5;
  cig.cigar_type = 'I';
  cigar_list.cigars.push_back(cig);
  cig.num = 90;
  cig.cigar_type = 'M';
  cigar_list.cigars.push_back(cig);
  cigar_list.ResetString();
  cig.num = 2;
  cig.cigar_type = 'I';
  cigar_list.cigars.push_back(cig);
  CPPUNIT_ASSERT_MESSAGE("ExtractCigar failed", !(ExtractCigar(cigar_list, 0, 80, 200, &str_cigar_list))); // boundary with indel at end
  CPPUNIT_ASSERT_MESSAGE("ExtractCigar failed", !(ExtractCigar(cigar_list, 0, 0, 50, &str_cigar_list))); // boundary with indel at end

  // Test invalid entries
  CPPUNIT_ASSERT_MESSAGE("ExtractCigar failed", !(ExtractCigar(cigar_list, 0, -1, 200, &str_cigar_list))); // invalid region start
  CPPUNIT_ASSERT_MESSAGE("ExtractCigar failed", !(ExtractCigar(cigar_list, 0, 100, 250, &str_cigar_list))); // invalid region end
  CPPUNIT_ASSERT_MESSAGE("ExtractCigar failed", !(ExtractCigar(cigar_list, -1, 100, 250, &str_cigar_list))); // invalid cigar start
  CPPUNIT_ASSERT_MESSAGE("ExtractCigar failed", !(ExtractCigar(cigar_list, 0, 100, 50, &str_cigar_list))); // start < end

  // Additional test cases 222M4I70M
  cigar_list.cigars.clear();
  cig.num = 222;
  cig.cigar_type = 'M';
  cigar_list.cigars.push_back(cig);
  cig.num = 4;
  cig.cigar_type = 'I';
  cigar_list.cigars.push_back(cig);
  cig.num = 70;
  cig.cigar_type = 'M';
  cigar_list.cigars.push_back(cig);
  expected_cigar_string = "222M4I70M";
  ExtractCigar(cigar_list, 18392754, 18392976, 18393035, &str_cigar_list);
  CPPUNIT_ASSERT_MESSAGE("ExtractCigar failed", str_cigar_list.cigar_string == expected_cigar_string);
}
