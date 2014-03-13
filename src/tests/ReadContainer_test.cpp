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

#include "src/tests/ReadContainer_test.h"
#include "src/runtime_parameters.h"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(ReadContainerTest);

using namespace std;
using BamTools::BamAlignment;

void ReadContainerTest::setUp() {
  vector<string> filenames;

  // This environment variable is defined in './src/Makefile.am',
  // Will be set during autotools' "make check" process.
  char* test_dir_env = getenv("LOBSTR_TEST_DIR");
  if (test_dir_env!=NULL)
    test_dir = test_dir_env;
  else
    test_dir = "../tests" ; // use as fall back
  string file = test_dir + "/test.aligned.sorted.bam" ;

  filenames.push_back(file);
  _read_container = new ReadContainer(filenames);

  // Don't check whether reads span STR
  min_border = -100000;
}

void ReadContainerTest::tearDown() {
  delete _read_container;
}

void ReadContainerTest::test_AddReadsFromFile() {


  map<pair<string,int>, string> ref_ext_nucleotides;
  ReferenceSTR ref_str;
  // When region is in the file
  ref_str.chrom = "chr22";
  ref_str.start = 38110651;
  ref_str.stop  = 38110651;
  _read_container->ClearReads();
  _read_container->AddReadsFromFile(ref_str, ref_ext_nucleotides);
  CPPUNIT_ASSERT_EQUAL(static_cast<int>(_read_container->aligned_str_map_.size()), 1);
  // When region not in file
  _read_container->ClearReads();
  ref_str.start = 0;
  ref_str.stop = 0;
  _read_container->AddReadsFromFile(ref_str, ref_ext_nucleotides);
  CPPUNIT_ASSERT(_read_container->aligned_str_map_.empty());
}

void ReadContainerTest::test_ParseRead() {
  map<pair<string,int>, string> ref_ext_nucleotides;
  BamAlignment aln;
  std::string rg = "test";
  std::string repseq = "AC";
  AlignedRead aligned_read;
  float copynum = 10;
  include_flank = false;
  // Test valid allele length
  aln.Name = "test";
  aln.QueryBases = "NNNNN";
  aln.Qualities = "NNNNN";
  aln.SetIsReverseStrand(true);
  aln.Position = 0;
  aln.SetIsSecondMate(false);
  aln.AddTag("RG","Z",rg);
  aln.AddTag("XS","i",0);
  aln.AddTag("XE","i",20);
  aln.AddTag("XR", "Z", repseq);
  aln.AddTag("XC", "f", copynum);
  aln.RefID = 0;
  // No XD
  aln.RemoveTag("XD");
  CPPUNIT_ASSERT(!(_read_container->ParseRead(aln, &aligned_read, ref_ext_nucleotides))); 
  // Test more valid allele lengths
  aln.AddTag("XD","i",20);
  CPPUNIT_ASSERT(_read_container->ParseRead(aln, &aligned_read, ref_ext_nucleotides));
  aln.RemoveTag("XD");
  aln.AddTag("XD","i",-19);
  CPPUNIT_ASSERT(_read_container->ParseRead(aln, &aligned_read, ref_ext_nucleotides));
  // Test invalid allele length
  aln.RemoveTag("XD");
  aln.AddTag("XD","i",-31);
  CPPUNIT_ASSERT(!(_read_container->ParseRead(aln, &aligned_read, ref_ext_nucleotides)));
  // Exceed max_diff_ref
  aln.RemoveTag("XD");
  aln.AddTag("XD","i",100);
  CPPUNIT_ASSERT(!(_read_container->ParseRead(aln, &aligned_read, ref_ext_nucleotides)));
  aln.RemoveTag("XD");
  aln.AddTag("XD","i",0);
  // Exceed max mate dist
  aln.AddTag("XM","i",1000000);
  CPPUNIT_ASSERT(!(_read_container->ParseRead(aln, &aligned_read, ref_ext_nucleotides)));
  aln.RemoveTag("XM");
  // Exceed max mapq
  aln.AddTag("XQ","i",10000);
  CPPUNIT_ASSERT(!(_read_container->ParseRead(aln, &aligned_read, ref_ext_nucleotides)));
  aln.RemoveTag("XQ");
  // Read is mate
  aln.SetIsSecondMate(true);
  CPPUNIT_ASSERT(!(_read_container->ParseRead(aln, &aligned_read, ref_ext_nucleotides)));
  aln.SetIsSecondMate(false);
  // Read is partial
  aln.AddTag("XP","i",1);
  CPPUNIT_ASSERT(!(_read_container->ParseRead(aln, &aligned_read, ref_ext_nucleotides)));
  aln.RemoveTag("XP");
  // Non-unit
  unit = true;
  aln.RemoveTag("XD");
  aln.AddTag("XD","i",5);
  CPPUNIT_ASSERT(!(_read_container->ParseRead(aln, &aligned_read, ref_ext_nucleotides)));
  unit = false;
}

void ReadContainerTest::test_GetReadsAtCoord() {
  _read_container->ClearReads();
  list<AlignedRead> results;
  pair<string, int> coord1 ("chr1", 0);
  _read_container->GetReadsAtCoord(coord1, &results);
  CPPUNIT_ASSERT_EQUAL(static_cast<int>(results.size()), 0);
  list<AlignedRead> list1;
  AlignedRead r1;
  AlignedRead r2;
  list1.push_back(r1);
  list1.push_back(r2);
  _read_container->aligned_str_map_.insert(pair<pair<string, int>, list<AlignedRead> >
					   (coord1, list1));
  _read_container->GetReadsAtCoord(coord1, &results);
  CPPUNIT_ASSERT_EQUAL(static_cast<int>(results.size()), 2);
  pair<string, int> coord2("chrX", 1);
  _read_container->GetReadsAtCoord(coord2, &results);
  CPPUNIT_ASSERT_EQUAL(static_cast<int>(results.size()), 0);
}

void ReadContainerTest::test_GetBamTags() {
  // TODO
  BamAlignment aln;
  string str_result;
  float float_result;
  string rg = "test";
  float xc = 12.0;
  aln.AddTag("RG","Z",rg);
  aln.AddTag("XD","i",0);
  aln.AddTag("XC","f",xc);
  // GetStringBamTag
  CPPUNIT_ASSERT(_read_container->GetStringBamTag(aln, "RG", &str_result));
  CPPUNIT_ASSERT(str_result == "test");
  CPPUNIT_ASSERT(!_read_container->GetStringBamTag(aln, "XX", &str_result));
  // GetFloatBamTag
  CPPUNIT_ASSERT(_read_container->GetFloatBamTag(aln, "XC", &float_result));
  CPPUNIT_ASSERT(float_result == 12.0);
  CPPUNIT_ASSERT(!_read_container->GetFloatBamTag(aln, "XX", &float_result));
  // GetIntBamTag
  int int_result;
  int d = 0;
  aln.AddTag("XA","i",d);
  CPPUNIT_ASSERT(_read_container->GetIntBamTag(aln, "XA", &int_result));
  CPPUNIT_ASSERT_EQUAL(int_result, 0);
  aln.RemoveTag("XA");
  int8_t d8 = 10;
  aln.AddTag("XA","c", d8);
  CPPUNIT_ASSERT(_read_container->GetIntBamTag(aln, "XA", &int_result));
  CPPUNIT_ASSERT_EQUAL(int_result, 10);
  aln.RemoveTag("XA");
  aln.AddTag("XA","c",rg);
  CPPUNIT_ASSERT(!_read_container->GetFloatBamTag(aln, "XA", &float_result));
  aln.RemoveTag("XA");
  uint8_t ud8 = 9;
  aln.AddTag("XA","C", ud8);
  CPPUNIT_ASSERT(_read_container->GetIntBamTag(aln, "XA", &int_result));
  CPPUNIT_ASSERT_EQUAL(int_result, 9);
  aln.RemoveTag("XA");
  aln.AddTag("XA","C",rg);
  CPPUNIT_ASSERT(!_read_container->GetFloatBamTag(aln, "XA", &float_result));
  aln.RemoveTag("XA");
  int16_t d16 = 8;
  aln.AddTag("XA","s",d16);
  CPPUNIT_ASSERT(_read_container->GetIntBamTag(aln, "XA", &int_result));
  CPPUNIT_ASSERT_EQUAL(int_result, 8);
  aln.RemoveTag("XA");
  aln.AddTag("XA","s",rg);
  CPPUNIT_ASSERT(!_read_container->GetFloatBamTag(aln, "XA", &float_result));
  aln.RemoveTag("XA");
  uint16_t ud16 = 7;
  aln.AddTag("XA","S",ud16);
  CPPUNIT_ASSERT(_read_container->GetIntBamTag(aln, "XA", &int_result));
  CPPUNIT_ASSERT_EQUAL(int_result, 7);
  aln.RemoveTag("XA");
  aln.AddTag("XA","S",rg);
  CPPUNIT_ASSERT(!_read_container->GetFloatBamTag(aln, "XA", &float_result));
  aln.RemoveTag("XA");
  uint32_t ud32 = 6;
  aln.AddTag("XA","I",ud32);
  CPPUNIT_ASSERT(_read_container->GetIntBamTag(aln, "XA", &int_result));
  CPPUNIT_ASSERT_EQUAL(int_result, 6);
  aln.RemoveTag("XA");
  aln.AddTag("XA","I",rg);
  CPPUNIT_ASSERT(!_read_container->GetFloatBamTag(aln, "XA", &float_result));
  aln.RemoveTag("XA");
}

void ReadContainerTest::test_GetSTRAllele() {
  CIGAR_LIST cigar_list;
  CIGAR cig1;
  cig1.num = 50;
  cig1.cigar_type = 'M';
  CIGAR cig2;
  cig2.num = 20;
  cig2.cigar_type = 'D';
  CIGAR cig3;
  cig3.num = 10;
  cig3.cigar_type = 'I';
  cigar_list.cigars.push_back(cig1);
  cigar_list.cigars.push_back(cig2);
  cigar_list.cigars.push_back(cig3);
  CPPUNIT_ASSERT_EQUAL(_read_container->GetSTRAllele(cigar_list), -10);
}
