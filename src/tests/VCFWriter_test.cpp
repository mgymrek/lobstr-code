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

#include "src/tests/VCFWriter_test.h"

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(VCFWriterTest);

void VCFWriterTest::setUp() {
  _vcf_writer = new VCFWriter("/dev/null");
}

void VCFWriterTest::tearDown() {
  delete _vcf_writer;
}

void VCFWriterTest::test_GetSTRVar() {
  // Case 1: reference
  std::string refseq = "TTTATTTATTTATTTATTTATTTATTTATTT";
  std::string ref_repseq = "TTTA";
  int allele = 0;
  CPPUNIT_ASSERT_EQUAL(refseq, _vcf_writer->GetSTRVar(refseq, ref_repseq, allele));

  // Case 2: < reference
  allele = -4;
  std::string check = "TTTATTTATTTATTTATTTATTTATTT";
  CPPUNIT_ASSERT_EQUAL(check, _vcf_writer->GetSTRVar(refseq, ref_repseq, allele));
  allele = -3;
  check = "TTTATTTATTTATTTATTTATTTATTTA";
  CPPUNIT_ASSERT_EQUAL(check, _vcf_writer->GetSTRVar(refseq, ref_repseq, allele));
  allele = -2;
  check = "TTTATTTATTTATTTATTTATTTATTTAT";
  CPPUNIT_ASSERT_EQUAL(check, _vcf_writer->GetSTRVar(refseq, ref_repseq, allele));
  allele = -1;
  check = "TTTATTTATTTATTTATTTATTTATTTATT";
  CPPUNIT_ASSERT_EQUAL(check, _vcf_writer->GetSTRVar(refseq, ref_repseq, allele));

  // Case 3: > reference
  allele = 1;
  check = "TTTATTTATTTATTTATTTATTTATTTATTTA";
  CPPUNIT_ASSERT_EQUAL(check, _vcf_writer->GetSTRVar(refseq, ref_repseq, allele));
  allele = 2;
  check = "TTTATTTATTTATTTATTTATTTATTTATTTAT";
  CPPUNIT_ASSERT_EQUAL(check, _vcf_writer->GetSTRVar(refseq, ref_repseq, allele));
  allele = 3;
  check = "TTTATTTATTTATTTATTTATTTATTTATTTATT";
  CPPUNIT_ASSERT_EQUAL(check, _vcf_writer->GetSTRVar(refseq, ref_repseq, allele));
  allele = 4;
  check = "TTTATTTATTTATTTATTTATTTATTTATTTATTT";
  CPPUNIT_ASSERT_EQUAL(check, _vcf_writer->GetSTRVar(refseq, ref_repseq, allele));
  allele = 5;
  check = "TTTATTTATTTATTTATTTATTTATTTATTTATTTA";
  CPPUNIT_ASSERT_EQUAL(check, _vcf_writer->GetSTRVar(refseq, ref_repseq, allele));
}
