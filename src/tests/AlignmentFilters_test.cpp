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
#include <cstdlib>

#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#include "src/tests/AlignmentFilters_test.h"
#include "src/runtime_parameters.h"

#include "common.h"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(AlignmentFiltersTest);

using namespace std;
using BamTools::BamAlignment;

vector<BamTools::CigarOp> GetCigarOps(string cigar_string){
  vector<BamTools::CigarOp> cigar_ops;
  unsigned int index = 0;
  while (index < cigar_string.size()){
    unsigned int num_start = index;
    while (cigar_string[index] >= '0' && cigar_string[index] <= '9' && index < cigar_string.size())
      index++;

    if (index == cigar_string.size() || index == num_start) PrintMessageDieOnError("Improperly formatted CIGAR string: " + cigar_string, ERROR);
    long num  = atol(cigar_string.substr(num_start, index-num_start).c_str());
    if (num <= 0) PrintMessageDieOnError("Improperly formatted CIGAR string: " + cigar_string, ERROR);
    char type = cigar_string[index];
    cigar_ops.push_back(BamTools::CigarOp(type, num));
    index++;
  }
  return cigar_ops;
}


void AlignmentFiltersTest::setUp() {
  
}

void AlignmentFiltersTest::tearDown() {
  
}

void AlignmentFiltersTest::test_GetDistToIndel(){
  AlignedRead aln;
  vector<BamTools::CigarOp> cigar_ops; 
  pair<int,int> dists;
  
  aln.cigar_ops = GetCigarOps("1M");
  dists         = AlignmentFilters::GetEndDistToIndel(&aln);
  CPPUNIT_ASSERT(dists.first == -1 && dists.second == -1);

  aln.cigar_ops = GetCigarOps("1D");
  dists         = AlignmentFilters::GetEndDistToIndel(&aln);
  CPPUNIT_ASSERT(dists.first == 0 && dists.second == 0);

  aln.cigar_ops = GetCigarOps("5H3S1M2I");
  dists         = AlignmentFilters::GetEndDistToIndel(&aln);
  CPPUNIT_ASSERT(dists.first == 1 && dists.second == 0);

  aln.cigar_ops = GetCigarOps("5H3S1M5D33M7S");
  dists         = AlignmentFilters::GetEndDistToIndel(&aln);
  CPPUNIT_ASSERT(dists.first == 1 && dists.second == 33);
}

void AlignmentFiltersTest::test_GetNumEndMatches(){
  AlignedRead aln;
  string ref_seq;
  int    ref_start;
  pair<int,int> num_matches;

  aln.cigar_ops   = GetCigarOps("5M1I3M");
  aln.read_start  = 100;
  aln.nucleotides = "ACACATCAC";
  ref_seq         = "ACACACACACACACACAC";
  ref_start       = 100;
  num_matches     = AlignmentFilters::GetNumEndMatches(&aln, ref_seq, ref_start);
  CPPUNIT_ASSERT(num_matches.first == 5 && num_matches.second == 3);

  aln.nucleotides = "ACACTTCTC";
  num_matches     = AlignmentFilters::GetNumEndMatches(&aln, ref_seq, ref_start);
  CPPUNIT_ASSERT(num_matches.first == 4 && num_matches.second == 1);
    
  aln.nucleotides = "ACACATCAC";
  ref_start       = 86;
  num_matches     = AlignmentFilters::GetNumEndMatches(&aln, ref_seq, ref_start);
  CPPUNIT_ASSERT(num_matches.first == -1 && num_matches.second == -1);

  ref_start   = 102;
  num_matches = AlignmentFilters::GetNumEndMatches(&aln, ref_seq, ref_start);
  CPPUNIT_ASSERT(num_matches.first == -1 && num_matches.second == -1);

  aln.cigar_ops   = GetCigarOps("5H5M1I3M3H");
  aln.read_start  = 100;
  aln.nucleotides = "ACACATCAC";
  ref_seq         = "ACACACACACACACACAC";
  ref_start       = 100;
  num_matches     = AlignmentFilters::GetNumEndMatches(&aln, ref_seq, ref_start);
  CPPUNIT_ASSERT(num_matches.first == 5 && num_matches.second == 3);

  aln.cigar_ops   = GetCigarOps("5H3S5M1I3M");
  aln.read_start  = 100;
  aln.nucleotides = "NNNACACATCAC";
  ref_seq         = "ACACACACACACACACAC";
  ref_start       = 100;
  num_matches     = AlignmentFilters::GetNumEndMatches(&aln, ref_seq, ref_start);
  CPPUNIT_ASSERT(num_matches.first == 5 && num_matches.second == 3);
}


