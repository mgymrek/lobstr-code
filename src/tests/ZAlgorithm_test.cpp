
#include <string>
#include <vector>

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
  for (int i = 0; i < NUM_TRIALS; i++){
    s1 = DNATools::RandDNA(LENGTH);
    s2 = DNATools::RandDNA(LENGTH);
    ZAlgorithm::GetPrefixMatchCounts(s1, s2, num_matches);
    for (unsigned int j = 0; j < s2.size(); j++){
      bool match  = (s1.substr(0, num_matches[j]).compare(s2.substr(j, num_matches[j])) == 0);
      bool longer = (s1.size() > num_matches[j] && s2.size() > (j + num_matches[j]) && s2[j+num_matches[j]] == s1[num_matches[j]]);
      CPPUNIT_ASSERT_MESSAGE("ZAlgorithm prefix error", !(!match || (match && longer)));
    }
  }
}

void ZAlgorithmTest::test_Suffix(){
  vector<int> num_matches;
  string s1, s2;
  for (int i = 0; i < NUM_TRIALS; i++){
    s1 = DNATools::RandDNA(LENGTH);
    s2 = DNATools::RandDNA(LENGTH);
    ZAlgorithm::GetSuffixMatchCounts(s1, s2, num_matches);
    for (int j = s2.size()-1; j >= 0; j--){
      string sub_1 = s1.substr(s1.size()-num_matches[j], num_matches[j]);
      string sub_2 = s2.substr(j-num_matches[j]+1, num_matches[j]);
      bool match   = (sub_1.compare(sub_2) == 0);
      bool longer  = ((j-num_matches[j]) > 0 && (s1.size()-num_matches[j]) > 0 && s1[s1.size()-num_matches[j]-1] == s2[j-num_matches[j]]);
      CPPUNIT_ASSERT_MESSAGE("ZAlgorithm prefix error", !(!match || (match && longer)));
    }
  }
}




