// Test functions in GST.h

#ifndef GENOTYPER_TEST_H_
#define GENOTYPER_TEST_H_

#include <cppunit/extensions/HelperMacros.h>
#include "Genotyper.h"

class GenotyperTest : public CppUnit::TestFixture { 
 public:
  static CppUnit::Test* suite();

  void setUp();

  void tearDown();

 protected:
  void testGetGenotype();
};

#endif /* GENOTYPER_TEST_H_ */
