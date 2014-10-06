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
#include "src/tests/AlignmentUtils_test.h"
#include "src/cigar.h"
#include "src/MSReadRecord.h"
#include "src/ReadPair.h"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(AlignmentUtilsTest);

void AlignmentUtilsTest::setUp() {}

void AlignmentUtilsTest::tearDown() {}

void AlignmentUtilsTest::test_StitchReads() {
  ReadPair read_pair;
  read_pair.aligned_read_num = 0;
  ALIGNMENT left_alignment;
  ALIGNMENT right_alignment;
  MSReadRecord read1;
  // Case 1: reads stitch nicely
  read1.orig_nucleotides = "AACCTACAGAAAATATTTCAGTATTTTGGTACATGAAGATCCATGTTGAGAGAGAAAGAAACAGAGAGTGAGAGAGAGACAGAGAGAGAGAGAGAGATTTG";
  read1.orig_qual = "ffffefffffee``ceeedeeeeeeeeeeaeeeeeee`eeeeeeeedeeedadeeeeeeedde`deeededaaedededadeaecdad`d_dadaba\\Y^`";
  MSReadRecord read2;
  read2.orig_nucleotides = "NAAGNACATTGGAAGTTTCTNTTCCNAATCTCTCTCTCTCTCTCTGTCTCTCTCTCACTCTCTGTTTCTTTCTCTCTCAACATGGATCTTCATGTACCAAA";
  read2.orig_qual = "BQQIBIIGJJ_``__BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB";
  read1.right_flank_index_from_end = 0;
  read1.left_flank_index_from_start = 0;
  read2.right_flank_index_from_end = 0;
  read2.left_flank_index_from_start = 0;
  left_alignment.left = true;
  read_pair.reads.push_back(read1);
  read_pair.reads.push_back(read2);
  CPPUNIT_ASSERT_MESSAGE("Stitching failed", AlignmentUtils::StitchReads(&read_pair, &left_alignment, &right_alignment));
  CPPUNIT_ASSERT_MESSAGE("Incorrect stitching returned - nucs", read_pair.reads.at(0).nucleotides == "AACCTACAGAAAATATTTCAGTATTTTGGTACATGAAGATCCATGTTGAGAGAGAAAGAAACAGAGAGTGAGAGAGAGACAGAGAGAGAGAGAGAGATTTGGAANAGAAACTTCCAATGTNCTTN");
  CPPUNIT_ASSERT_MESSAGE("Incorrect stitching returned - qual", read_pair.reads.at(0).quality_scores == "ffffefffffee``ceeedeeeeeeeeeeaeeeeeee`eeeeeeeedeeedadeeeeeeedde`deeededaaedededadeaecdad`d_dadaba\\Y^`BBBBBBBBB__``_JJGIIBIQQB");
  read_pair.reads.at(0).orig_nucleotides = "ACGTAATTAATAATAATAATAATAATAATAATAATAATAATAATAAAGTAGCCAGGTATGGAGGCACAGGTCTGTAGTACCAGCTG";
  read_pair.reads.at(1).orig_nucleotides = "TGAGCTCAAGTTGTCCTTCTGCTTCAGCTTCCCAAGTAGCTGGGACTACAGACCTGTGCCTCCATACCTGGCTACTTTATTATTATTATTATTATTATTAT";
  read_pair.reads.at(0).orig_qual = "26:=67;>:;?>:?><>>=<A=>9<<>=?>><;2996;98<<>:B>?<>A>:44?;=?;<;=;8585:8:73578:3461222010";
  read_pair.reads.at(1).orig_qual = ".4339638<4=7;==;@>=@?<B@?@><?=::;;>::?;8=<?A>:==:=?:=7<:<>96767997456=:55516886:7474685356174443/5031";
  CPPUNIT_ASSERT_MESSAGE("Stitching failed", AlignmentUtils::StitchReads(&read_pair, &left_alignment, &right_alignment));
  CPPUNIT_ASSERT_MESSAGE("Incorrect stitching returned - nucs", read_pair.reads.at(0).nucleotides == "ACGTAATTAATAATAATAATAATAATAATAATAATAATAATAATAAAGTAGCCAGGTATGGAGGCACAGGTCTGTAGTCCCAGCTACTTGGGAAGCTGAAGCAGAAGGACAACTTGAGCTCA");

  // Case 1.5 reads stitch nicely from other direction
  read_pair.reads.at(0).orig_nucleotides = "TGTATTTCATGTGTACATTCGTATCTATCTATCTATCTATCTATCTATCCATCTATCTATCTATTCCCCACAGTGAAAATAATCTACAGGATAGGTAAATA";
  read_pair.reads.at(1).orig_nucleotides = "ATTTACCTATCCTGTAGATTATTTTCACTGTGGGGAATAGATAGATAGATGGATAGATAGATAGATAGATAGATAGATACGAATGTACACATGAAATACAA";
read_pair.reads.at(0).orig_qual = "TGTATTTCATGTGTACATTCGTATCTATCTATCTATCTATCTATCTATCCATCTATCTATCTATTCCCCACAGTGAAAATAATCTACAGGATAGGTAAATA";
  read_pair.reads.at(1).orig_qual = "ATTTACCTATCCTGTAGATTATTTTCACTGTGGGGAATAGATAGATAGATGGATAGATAGATAGATAGATAGATAGATACGAATGTACACATGAAATACAA";
  CPPUNIT_ASSERT_MESSAGE("Stitching failed", AlignmentUtils::StitchReads(&read_pair, &left_alignment, &right_alignment));
  CPPUNIT_ASSERT_MESSAGE("Incorrect stitching returned - nucs", read_pair.reads.at(0).nucleotides == "TTGTATTTCATGTGTACATTCGTATCTATCTATCTATCTATCTATCTATCCATCTATCTATCTATTCCCCACAGTGAAAATAATCTACAGGATAGGTAAATA");

  // Case 2: too repetitive to stitch
  read_pair.reads.at(0).orig_nucleotides = "AACCTACAGAAAATATTTCAGTATTTTGGTACATGAAGATCCATGTTGCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCT";
  read_pair.reads.at(0).orig_qual = "AACCTACAGAAAATATTTCAGTATTTTGGTACATGAAGATCCATGTTGAGAGAGAAAGAAACAGAGAGTGAAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG";
  read_pair.reads.at(1).orig_nucleotides = "GTTTCTTTCTCTCTCAACATGGATCTTCATGTACCAAAATACTGAAATATTTTCTGTAGGTTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCT";
  read_pair.reads.at(1).orig_qual = "GTTTCTTTCTCTCTCAACATGGATCTTCATGTACCAAAATACTGAAATATTTTCTGTAGGTTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCT";
  CPPUNIT_ASSERT_MESSAGE("This stitch should be too repetitive", !(AlignmentUtils::StitchReads(&read_pair, &left_alignment, &right_alignment)));

  // Case 3: reads don't stitch
  read_pair.reads.at(0).orig_nucleotides = "AACCTACAGAAAATATTTCAGTATTTTGGTACATGAAGATCCATGTTGAGAGAGAAAGAAACAGAGAGTGAAGAGAGAGAGAGAGAGAGAGAGAGAGAGAG";
  read_pair.reads.at(1).orig_nucleotides = "AGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAACCTACAGAAAATATTTCAGTATTTTGGTACATGAAGATCCATGTTGAGAGAGAAAGAAAC";
  CPPUNIT_ASSERT_MESSAGE("This should not stitch", !(AlignmentUtils::StitchReads(&read_pair, &left_alignment, &right_alignment)));
}

void AlignmentUtilsTest::test_GetMapq() {
  int edit;
  // Case 1: no mismatches
  std::string aligned = "ACACGTACGTTATCGATCGAT";
  std::string ref = "ACACGTACGTTATCGATCGAT";
  std::string qual = "eeeeedeeedadeeeeeeedd";
  CPPUNIT_ASSERT_EQUAL(0, AlignmentUtils::GetMapq(aligned, ref, qual, &edit));
  CPPUNIT_ASSERT_EQUAL(edit, 0);

  // Case 2: mismatches
  aligned = "TCACGTACGTTATCGATCGAT";
  CPPUNIT_ASSERT_EQUAL(static_cast<int>('e')-33, AlignmentUtils::GetMapq(aligned, ref, qual, &edit));
  CPPUNIT_ASSERT_EQUAL(edit, 1);

  // Case 3: gaps, no mismatches
  aligned = "AC-ACGTACGTTATCGATCGAT";
  ref = "ACTACGTACGTTATCGATCGAT";
  CPPUNIT_ASSERT_EQUAL(0, AlignmentUtils::GetMapq(aligned, ref, qual, &edit));
  CPPUNIT_ASSERT_EQUAL(edit, 0);

  // Case 4: gaps and mismatches
  aligned = "AC-ACGTACGTTATCGATCGAC";
  CPPUNIT_ASSERT_EQUAL(static_cast<int>('d')-33, AlignmentUtils::GetMapq(aligned, ref, qual, &edit));
  CPPUNIT_ASSERT_EQUAL(edit, 1);
}

void AlignmentUtilsTest::test_GetSTRAllele() {
  // Set up read record and test cigar
  MSReadRecord testRead;
  testRead.msStart = 1000;
  testRead.msEnd = 1050;
  testRead.read_start = 950;
  testRead.nucleotides = "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN";
  CIGAR_LIST testCigar;

  // Case 1: All match reference
  testCigar.cigars.clear();
  CIGAR match_cigar;
  match_cigar.cigar_type='M';
  match_cigar.num = 150;
  testCigar.cigars.push_back(match_cigar);
  testCigar.ResetString();
  CPPUNIT_ASSERT_MESSAGE("Incorrect CIGAR string", "150M" == testCigar.cigar_string);
  CPPUNIT_ASSERT_MESSAGE("GetSTRAllele failed", AlignmentUtils::GetSTRAllele(&testRead, testCigar));
  CPPUNIT_ASSERT_EQUAL(0, testRead.diffFromRef);


  // Case 2: Del part of STR
  // 2a: del at beginning of STR
  testCigar.cigars.clear();
  CIGAR l_cigar, str_cigar, r_cigar;
  l_cigar.cigar_type = 'M';
  l_cigar.num = 50;
  str_cigar.cigar_type = 'D';
  str_cigar.num = 4;
  r_cigar.cigar_type = 'M';
  r_cigar.num = 100;
  testCigar.cigars.push_back(l_cigar);
  testCigar.cigars.push_back(str_cigar);
  testCigar.cigars.push_back(r_cigar);
  testCigar.ResetString();
  CPPUNIT_ASSERT_MESSAGE("Incorrect CIGAR string", "50M4D100M" == testCigar.cigar_string);
  CPPUNIT_ASSERT_MESSAGE("GetSTRAllele failed", AlignmentUtils::GetSTRAllele(&testRead, testCigar));
  CPPUNIT_ASSERT_EQUAL(-4, testRead.diffFromRef);

  // 2b: del at end of STR
  testCigar.cigars.clear();
  l_cigar.cigar_type = 'M';
  l_cigar.num = 100;
  str_cigar.cigar_type = 'D';
  str_cigar.num = 4;
  r_cigar.cigar_type = 'M';
  r_cigar.num = 50;
  testCigar.cigars.push_back(l_cigar);
  testCigar.cigars.push_back(str_cigar);
  testCigar.cigars.push_back(r_cigar);
  testCigar.ResetString();
  CPPUNIT_ASSERT_MESSAGE("Incorrect CIGAR string", "100M4D50M" == testCigar.cigar_string);
  CPPUNIT_ASSERT_MESSAGE("GetSTRAllele failed", AlignmentUtils::GetSTRAllele(&testRead, testCigar));
  CPPUNIT_ASSERT_EQUAL(-4, testRead.diffFromRef);

  // 2c: del in middle of STR
  testCigar.cigars.clear();
  l_cigar.cigar_type = 'M';
  l_cigar.num = 75;
  str_cigar.cigar_type = 'D';
  str_cigar.num = 4;
  r_cigar.cigar_type = 'M';
  r_cigar.num = 75;
  testCigar.cigars.push_back(l_cigar);
  testCigar.cigars.push_back(str_cigar);
  testCigar.cigars.push_back(r_cigar);
  testCigar.ResetString();
  CPPUNIT_ASSERT_MESSAGE("Incorrect CIGAR string", "75M4D75M" == testCigar.cigar_string);
  CPPUNIT_ASSERT_MESSAGE("GetSTRAllele failed", AlignmentUtils::GetSTRAllele(&testRead, testCigar));
  CPPUNIT_ASSERT_EQUAL(-4, testRead.diffFromRef);

  // Case 3: Ins part of STR
  // 3a: Ins at beginning of STR
  testCigar.cigars.clear();
  l_cigar.cigar_type = 'M';
  l_cigar.num = 50;
  str_cigar.cigar_type = 'I';
  str_cigar.num = 4;
  r_cigar.cigar_type = 'M';
  r_cigar.num = 96;
  testCigar.cigars.push_back(l_cigar);
  testCigar.cigars.push_back(str_cigar);
  testCigar.cigars.push_back(r_cigar);
  testCigar.ResetString();
  CPPUNIT_ASSERT_MESSAGE("Incorrect CIGAR string", "50M4I96M" == testCigar.cigar_string);
  CPPUNIT_ASSERT_MESSAGE("GetSTRAllele failed", AlignmentUtils::GetSTRAllele(&testRead, testCigar));
  CPPUNIT_ASSERT_EQUAL(4, testRead.diffFromRef);

  // 3b: Ins at end of STR
  testCigar.cigars.clear();
  l_cigar.cigar_type = 'M';
  l_cigar.num = 100;
  str_cigar.cigar_type = 'I';
  str_cigar.num = 4;
  r_cigar.cigar_type = 'M';
  r_cigar.num = 46;
  testCigar.cigars.push_back(l_cigar);
  testCigar.cigars.push_back(str_cigar);
  testCigar.cigars.push_back(r_cigar);
  testCigar.ResetString();
  CPPUNIT_ASSERT_MESSAGE("Incorrect CIGAR string", "100M4I46M" == testCigar.cigar_string);
  CPPUNIT_ASSERT_MESSAGE("GetSTRAllele failed", AlignmentUtils::GetSTRAllele(&testRead, testCigar));
  CPPUNIT_ASSERT_EQUAL(4, testRead.diffFromRef);

  // 3c: Ins in middle of STR
  testCigar.cigars.clear();
  l_cigar.cigar_type = 'M';
  l_cigar.num = 75;
  str_cigar.cigar_type = 'I';
  str_cigar.num = 4;
  r_cigar.cigar_type = 'M';
  r_cigar.num = 71;
  testCigar.cigars.push_back(l_cigar);
  testCigar.cigars.push_back(str_cigar);
  testCigar.cigars.push_back(r_cigar);
  testCigar.ResetString();
  CPPUNIT_ASSERT_MESSAGE("Incorrect CIGAR string", "75M4I71M" == testCigar.cigar_string);
  CPPUNIT_ASSERT_MESSAGE("GetSTRAllele failed", AlignmentUtils::GetSTRAllele(&testRead, testCigar));
  CPPUNIT_ASSERT_EQUAL(4, testRead.diffFromRef);


  // Case 4: Del part of L flank
  testCigar.cigars.clear();
  l_cigar.cigar_type = 'M';
  l_cigar.num = 25;
  str_cigar.cigar_type = 'D';
  str_cigar.num = 4;
  r_cigar.cigar_type = 'M';
  r_cigar.num = 125;
  testCigar.cigars.push_back(l_cigar);
  testCigar.cigars.push_back(str_cigar);
  testCigar.cigars.push_back(r_cigar);
  testCigar.ResetString();
  CPPUNIT_ASSERT_MESSAGE("Incorrect CIGAR string", "25M4D125M" == testCigar.cigar_string);
  CPPUNIT_ASSERT_MESSAGE("GetSTRAllele failed", AlignmentUtils::GetSTRAllele(&testRead, testCigar));
  CPPUNIT_ASSERT_EQUAL(0, testRead.diffFromRef);

  // Case 5: Del part of R flank
  testCigar.cigars.clear();
  l_cigar.cigar_type = 'M';
  l_cigar.num = 125;
  str_cigar.cigar_type = 'D';
  str_cigar.num = 4;
  r_cigar.cigar_type = 'M';
  r_cigar.num = 25;
  testCigar.cigars.push_back(l_cigar);
  testCigar.cigars.push_back(str_cigar);
  testCigar.cigars.push_back(r_cigar);
  testCigar.ResetString();
  CPPUNIT_ASSERT_MESSAGE("Incorrect CIGAR string", "125M4D25M" == testCigar.cigar_string);
  CPPUNIT_ASSERT_MESSAGE("GetSTRAllele failed", AlignmentUtils::GetSTRAllele(&testRead, testCigar));
  CPPUNIT_ASSERT_EQUAL(0, testRead.diffFromRef);

  // Case 6: Ins part of L flank
  testCigar.cigars.clear();
  l_cigar.cigar_type = 'M';
  l_cigar.num = 25;
  str_cigar.cigar_type = 'I';
  str_cigar.num = 4;
  r_cigar.cigar_type = 'M';
  r_cigar.num = 121;
  testCigar.cigars.push_back(l_cigar);
  testCigar.cigars.push_back(str_cigar);
  testCigar.cigars.push_back(r_cigar);
  testCigar.ResetString();
  CPPUNIT_ASSERT_MESSAGE("Incorrect CIGAR string", "25M4I121M" == testCigar.cigar_string);
  CPPUNIT_ASSERT_MESSAGE("GetSTRAllele failed", AlignmentUtils::GetSTRAllele(&testRead, testCigar));
  CPPUNIT_ASSERT_EQUAL(0, testRead.diffFromRef);

  // Case 7: Ins part of R flank
  testCigar.cigars.clear();
  l_cigar.cigar_type = 'M';
  l_cigar.num = 125;
  str_cigar.cigar_type = 'I';
  str_cigar.num = 4;
  r_cigar.cigar_type = 'M';
  r_cigar.num = 21;
  testCigar.cigars.push_back(l_cigar);
  testCigar.cigars.push_back(str_cigar);
  testCigar.cigars.push_back(r_cigar);
  testCigar.ResetString();
  CPPUNIT_ASSERT_MESSAGE("Incorrect CIGAR string", "125M4I21M" == testCigar.cigar_string);
  CPPUNIT_ASSERT_MESSAGE("GetSTRAllele failed", AlignmentUtils::GetSTRAllele(&testRead, testCigar));
  CPPUNIT_ASSERT_EQUAL(0, testRead.diffFromRef);

  // Case 8: Ins in all three
  testCigar.cigars.clear();
  CIGAR l1, l2, l3, str1, r1, r2, r3;
  l1.cigar_type = 'M';
  l1.num = 25;
  l2.cigar_type = 'I';
  l2.num = 2;
  l3.cigar_type = 'M';
  l3.num = 48;
  str1.cigar_type = 'I';
  str1.num = 2;
  r1.cigar_type = 'M';
  r1.num = 47;
  r2.cigar_type = 'I';
  r2.num = 3;
  r3.cigar_type = 'M';
  r3.num = 23;
  testCigar.cigars.push_back(l1);
  testCigar.cigars.push_back(l2);
  testCigar.cigars.push_back(l3);
  testCigar.cigars.push_back(str1);
  testCigar.cigars.push_back(r1);
  testCigar.cigars.push_back(r2);
  testCigar.cigars.push_back(r3);
  testCigar.ResetString();
  CPPUNIT_ASSERT_MESSAGE("Incorrect CIGAR string", "25M2I48M2I47M3I23M" == testCigar.cigar_string);

  // Case 9: Del in all three
  testCigar.cigars.clear();
  l1.cigar_type = 'M';
  l1.num = 25;
  l2.cigar_type = 'D';
  l2.num = 2;
  l3.cigar_type = 'M';
  l3.num = 50;
  str1.cigar_type = 'D';
  str1.num = 2;
  r1.cigar_type = 'M';
  r1.num = 50;
  r2.cigar_type = 'D';
  r2.num = 3;
  r3.cigar_type = 'M';
  r3.num = 25;
  testCigar.cigars.push_back(l1);
  testCigar.cigars.push_back(l2);
  testCigar.cigars.push_back(l3);
  testCigar.cigars.push_back(str1);
  testCigar.cigars.push_back(r1);
  testCigar.cigars.push_back(r2);
  testCigar.cigars.push_back(r3);
  testCigar.ResetString();
  CPPUNIT_ASSERT_MESSAGE("Incorrect CIGAR string", "25M2D50M2D50M3D25M" == testCigar.cigar_string);
}

