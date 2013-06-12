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

void CommonTest::test_getMSSeq() {
  // Case 1: simple
  std::string nucs = "ATATATATATATATATAT";
  int k = 2;
  std::string repeat;
  CPPUNIT_ASSERT_MESSAGE("getMSSeq failed", getMSSeq(nucs, k, &repeat));
  CPPUNIT_ASSERT_MESSAGE("wrong motif returned", repeat == "AT");

  // Case 2: mismatches
  nucs = "ATATATATATATGTATATAT";
  CPPUNIT_ASSERT_MESSAGE("getMSSeq failed", getMSSeq(nucs, k, &repeat));
  CPPUNIT_ASSERT_MESSAGE("wrong motif returned", repeat == "AT");
  
  // Case 3: indels
  nucs = "ATATATATATTATATATATAT";
  CPPUNIT_ASSERT_MESSAGE("getMSSeq failed", getMSSeq(nucs, k, &repeat));
  CPPUNIT_ASSERT_MESSAGE("wrong motif returned", repeat == "AT");
  
  // Case 4: both
  nucs = "CAGCAGCAGCAGCAGGCAGCAGCAT";
  k = 3;
  CPPUNIT_ASSERT_MESSAGE("getMSSeq failed", getMSSeq(nucs, k, &repeat));
  CPPUNIT_ASSERT_MESSAGE("wrong motif returned", repeat == "AGC");

  // Case 5: nucs too short
  nucs = "CCG";
  k = 6;
  CPPUNIT_ASSERT_MESSAGE("getMSSeq should have failed", !getMSSeq(nucs, k, &repeat));
}

void CommonTest::test_getMinPermutation() {
  CPPUNIT_ASSERT_MESSAGE("Wrong minimum permutation", getMinPermutation("A") == "A");
  CPPUNIT_ASSERT_MESSAGE("Wrong minimum permutation", getMinPermutation("G") == "G");
  CPPUNIT_ASSERT_MESSAGE("Wrong minimum permutation", getMinPermutation("TACG") == "ACGT");
  CPPUNIT_ASSERT_MESSAGE("Wrong minimum permutation", getMinPermutation("TAGTACTAT") == "ACTATTAGT");
  CPPUNIT_ASSERT_MESSAGE("Wrong minimum permutation", getMinPermutation("CGCTCCC") == "CCCCGCT");
}

void CommonTest::test_getCanonicalRepeat() {
  CPPUNIT_ASSERT_MESSAGE("Wrong canonical repeat", getCanonicalRepeat("AAAAAAA") == "A");
  CPPUNIT_ASSERT_MESSAGE("Wrong canonical repeat", getCanonicalRepeat("TTT") == "A");
  CPPUNIT_ASSERT_MESSAGE("Wrong canonical repeat", getCanonicalRepeat("ATCGATC") == "ATCATCG");
  CPPUNIT_ASSERT_MESSAGE("Wrong canonical repeat", getCanonicalRepeat("TATATA") == "AT");
  CPPUNIT_ASSERT_MESSAGE("Wrong canonical repeat", getCanonicalRepeat("AGTCAGTC") == "ACTG");
}


void CommonTest::test_getCanonicalMS() {
  std::string seq;
  getCanonicalMS("AAAAA", &seq);
  CPPUNIT_ASSERT_MESSAGE("Wrong canonical sequence", seq == "AAAAA");
  
  getCanonicalMS("C", &seq);
  CPPUNIT_ASSERT_MESSAGE("Wrong canonical sequence", seq == "C");

  getCanonicalMS("T", &seq);
  CPPUNIT_ASSERT_MESSAGE("Wrong canonical sequence", seq == "A");

  getCanonicalMS("GGG", &seq);
  CPPUNIT_ASSERT_MESSAGE("Wrong canonical sequence", seq == "CCC");

  getCanonicalMS("CGACG", &seq);
  CPPUNIT_ASSERT_MESSAGE("Wrong canonical sequence", seq == "ACGCG");

  getCanonicalMS("GCTGC", &seq);
  CPPUNIT_ASSERT_MESSAGE("Wrong canonical sequence", seq == "AGCGC");
}



void CommonTest::test_IsPerfectRepeat() {
  // Case 1: perfect repeat
  std::string nucs = "ATATATATATATATATATATAT";
  std::string repeat = "AT";
  CPPUNIT_ASSERT_MESSAGE("Should be perfect repeat", IsPerfectRepeat(nucs, repeat));

  // Case 2: obviously not perfect repeat
  nucs = "ACGATTGCGCGGCGG";
  repeat = "AT";
  CPPUNIT_ASSERT_MESSAGE("Not perfect", !IsPerfectRepeat(nucs, repeat));

  // Case 3: close but not perfect
  nucs = "ATATATAGATATATA";
  repeat = "AT";
  CPPUNIT_ASSERT_MESSAGE("Not perfect", !IsPerfectRepeat(nucs, repeat));
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
