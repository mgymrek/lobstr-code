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

#include <string>
#include <vector>

#include <stdlib.h>

#include "src/tests/DNATools.h"
#include "src/tests/ZAlgorithm_test.h"
#include "src/common.h"
#include "src/ZAlgorithm.h"


using namespace std;

void ZAlgorithmTest::setUp(){}
void ZAlgorithmTest::tearDown(){}

void ZAlgorithmTest::test_Prefix(){
  vector<int> num_matches;
  string s1, s2;
  size_t s2_start, s2_stop;
  for (int i = 0; i < NUM_TRIALS; i++){
    s1       = DNATools::RandDNA(1+rand()%(LENGTH-1));
    s2       = DNATools::RandDNA(1+rand()%(LENGTH-1));
    s2_start = rand()%s2.size();
    s2_stop  = s2_start + (rand()%(s2.size()-s2_start));
    ZAlgorithm::GetPrefixMatchCounts(s1, s2, s2_start, s2_stop, num_matches);
    for (size_t j = s2_start; j <= s2_stop; j++){
      bool match  = (s1.substr(0, (size_t)num_matches[j-s2_start]).compare(s2.substr(j, (size_t)num_matches[j-s2_start])) == 0);
      bool longer = (s1.size() > (size_t)num_matches[j-s2_start] && s2.size() > (j + (size_t)num_matches[j-s2_start]) && s2[j+(size_t)num_matches[j-s2_start]] == s1[(size_t)num_matches[j-s2_start]]);
      stringstream msg;
      msg << "ZAlgorithm suffix error: s1 " << s1 << " s2 " << s2 << " num_matches " << num_matches[j] << " j " << j;
      CPPUNIT_ASSERT_MESSAGE(msg.str(), !(!match || (match && longer)));
    }
  }
}

void ZAlgorithmTest::test_Suffix(){
  vector<int> num_matches;
  string s1, s2;
  size_t s2_start, s2_stop;
  for (int i = 0; i < NUM_TRIALS; i++){
    s1       = DNATools::RandDNA(1+rand()%(LENGTH-1));
    s2       = DNATools::RandDNA(1+rand()%(LENGTH-1));
    s2_start = rand()%s2.size();
    s2_stop  = s2_start + (rand()%(s2.size()-s2_start));
    ZAlgorithm::GetSuffixMatchCounts(s1, s2, s2_start, s2_stop, num_matches);
    for (size_t j = s2_start; j <= s2_stop; j++){
      string sub_1 = s1.substr(s1.size()-(size_t)num_matches[j-s2_start], (size_t)num_matches[j-s2_start]);
      string sub_2 = s2.substr(j-(size_t)num_matches[j-s2_start]+1, (size_t)num_matches[j-s2_start]);
      bool match   = (sub_1.compare(sub_2) == 0);
      bool longer  = ((j-(size_t)num_matches[j-s2_start]) > 0 && (s1.size()-(size_t)num_matches[j-s2_start]) > 0 && s1[s1.size()-(size_t)num_matches[j-s2_start]-1] == s2[j-(size_t)num_matches[j-s2_start]]);
      stringstream msg;
      msg << "ZAlgorithm suffix error: s1 " << s1 << " s2 " << s2 << " num_matches " << num_matches[j] << " j " << j;
      CPPUNIT_ASSERT_MESSAGE(msg.str(), !(!match || (match && longer)));
    }
  }
}




