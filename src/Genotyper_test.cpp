// Test functions in common.h

#include <string>

#include "Genotyper_test.h"

using namespace std;

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(GenotyperTest);

void GenotyperTest::setUp() {}

void GenotyperTest::tearDown() {}

CppUnit::Test* GenotyperTest::suite() { 
  CppUnit::TestSuite *suiteOfTests = 
    new CppUnit::TestSuite("GenotyperTest");
  suiteOfTests->
    addTest(new CppUnit::TestCaller<GenotyperTest>("testGetGenotype", 
						   &GenotyperTest::testGetGenotype) );
  return suiteOfTests;
}

void GenotyperTest::testGetGenotype() {
  CPPUNIT_FAIL("not implemented");
}
