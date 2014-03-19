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

#include <cppunit/ui/text/TestRunner.h>

#include "src/tests/AlignmentFilters_test.h"
#include "src/tests/AlignmentUtils_test.h"
#include "src/tests/common_test.h"
#include "src/tests/NWNoRefEndPenalty_test.h"
#include "src/tests/ReadContainer_test.h"
#include "src/tests/RemoveDuplicates_test.h"
#include "src/tests/VCFWriter_test.h"

int main( int argc, char **argv) {
  // Adds the test to the list of tests to run
  CppUnit::TextUi::TestRunner runner;
  runner.addTest(AlignmentFiltersTest::suite());
  runner.addTest(AlignmentUtilsTest::suite());
  runner.addTest(CommonTest::suite());
  runner.addTest(NWNoRefEndPenaltyTest::suite());
  runner.addTest(ReadContainerTest::suite());
  runner.addTest(RemoveDuplicatesTest::suite());
  runner.addTest(VCFWriterTest::suite());

  // Run the tests
  bool wasSucessful = runner.run();

  // Return error code 1 if the one of test failed.
  return wasSucessful ? 0 : 1;
}

