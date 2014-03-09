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
#include <cstdlib>

#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#include "src/AlignedRead.h"
#include "src/tests/RemoveDuplicates_test.h"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(RemoveDuplicatesTest);

using namespace std;
using RemoveDuplicates::GetScore;
using RemoveDuplicates::GetRepRead;
using RemoveDuplicates::RemovePCRDuplicates;

void RemoveDuplicatesTest::setUp() {}

void RemoveDuplicatesTest::tearDown() {}

void RemoveDuplicatesTest::test_GetScore() {
  string qualstring = "==="; // = is 61-33=28
  float score = 28;
  CPPUNIT_ASSERT_EQUAL(RemoveDuplicates::GetScore(qualstring), score);
  qualstring = "===AAA";
  score = 30;
  CPPUNIT_ASSERT_EQUAL(RemoveDuplicates::GetScore(qualstring), score);
}

void RemoveDuplicatesTest::test_GetRepRead() {
  list<AlignedRead> aligned_reads;
  AlignedRead read1, read2, rep_read;
  read1.qualities = "========";
  read2.qualities = "AAAAAAAA";
  aligned_reads.push_back(read1);
  aligned_reads.push_back(read2);
  RemoveDuplicates::GetRepRead(aligned_reads, &rep_read);
  CPPUNIT_ASSERT_EQUAL(rep_read.qualities, read2.qualities);
}

void RemoveDuplicatesTest::test_RemovePCRDuplicates() {
  list<AlignedRead> aligned_reads;
  AlignedRead read1a, read1b, read1c, read2a, read2b, read3;
  read1a.read_start = 10;
  read1a.nucleotides = "AAAAA";
  read1a.qualities = "=====";
  read1b.read_start = 10;
  read1b.nucleotides = "AAAAA";
  read1b.qualities = "=====";
  read1c.read_start = 10;
  read1c.nucleotides = "AAAAA";
  read1c.nucleotides = "AAAAA";
  read2a.read_start = 20;
  read2a.nucleotides = "AAAAA";
  read2a.qualities = "=====";
  read2b.read_start = 20;
  read2b.nucleotides = "AAAAA";
  read2b.qualities = "AAAAA";
  read3.read_start = 30;
  read3.nucleotides = "AAAAA";
  read3.qualities = "=====";
  aligned_reads.push_back(read1a);
  aligned_reads.push_back(read1b);
  aligned_reads.push_back(read1c);
  aligned_reads.push_back(read2a);
  aligned_reads.push_back(read2b);
  aligned_reads.push_back(read3);
  RemoveDuplicates::RemovePCRDuplicates(&aligned_reads);
  CPPUNIT_ASSERT_EQUAL(static_cast<int>(aligned_reads.size()), 3);
}
