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

#include "src/tests/ReadContainer_test.h"
#include "src/runtime_parameters.h"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(ReadContainerTest);

using namespace std;
using BamTools::BamAlignment;

void ReadContainerTest::setUp() {
  vector<string> filenames;
  filenames.push_back("../tests/test.aligned.sorted.bam"); // TODO how to get path
  _read_container = new ReadContainer(filenames);
}

void ReadContainerTest::tearDown() {
  delete _read_container;
}

void ReadContainerTest::test_ParseRead() {
  BamAlignment aln;
  std::string rg = "test";
  std::string repseq = "AC";
  AlignedRead aligned_read;
  float copynum = 25;
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
  aln.AddTag("XE","i",50);
  aln.AddTag("XD","i",0);
  aln.AddTag("XR", "Z", repseq);
  aln.AddTag("XC", "f", copynum);
  aln.RefID = 0;
  CPPUNIT_ASSERT(_read_container->ParseRead(aln, &aligned_read));
  // Test more valid allele lengths
  aln.RemoveTag("XD");
  aln.AddTag("XD","i",40);
  CPPUNIT_ASSERT(_read_container->ParseRead(aln, &aligned_read));
  aln.RemoveTag("XD");
  aln.AddTag("XD","i",-49);
  CPPUNIT_ASSERT(_read_container->ParseRead(aln, &aligned_read));
  // Test invalid allele length
  aln.RemoveTag("XD");
  aln.AddTag("XD","i",-51);
  CPPUNIT_ASSERT(!(_read_container->ParseRead(aln, &aligned_read)));
}
