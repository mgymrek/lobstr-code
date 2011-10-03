// Test suite for STR functions

#include <cppunit/ui/text/TestRunner.h>

#include "Genotyper_test.h"

int main( int argc, char **argv)
{
  // Adds the test to the list of tests to run
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(GenotyperTest::suite());
  // Run the tests.
  bool wasSucessful = runner.run();

  // Return error code 1 if the one of test failed.
  return wasSucessful ? 0 : 1;
}

