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
#include <sstream>

#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>

#include "src/common.h"
#include "src/tests/DNATools.h"
#include "src/tests/NWNoRefEndPenalty_test.h"

using namespace std;

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(NWNoRefEndPenaltyTest);


void NWNoRefEndPenaltyTest::setUp() {} 
void NWNoRefEndPenaltyTest::tearDown() {}

int NWNoRefEndPenaltyTest::GenAlignments(int num_trials, double mut_prob, double ins_prob, double del_prob){
  srand(RAND_SEED);
  int num_correct = 0;
  for (int i = 0; i < num_trials; i++){
    // Create random reference sequence
    string ref_seq = DNATools::RandDNA(REF_LEN);

    // Choose the read region
    int read_start  = rand()%(ref_seq.size()-READ_LEN+1);
    string read_seq = ref_seq.substr(read_start, READ_LEN);

    // Randomly mutate read region 
    stringstream mut_read_seq_ss;
    for (unsigned int i = PERFECT_FLANK; i < read_seq.size()-PERFECT_FLANK; i++){
      float prob = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
      if (prob <= mut_prob)
	mut_read_seq_ss << DNATools::RandDNA(1);
      else if (prob - mut_prob <= ins_prob){
	int ins_size = rand()%(MAX_INS-1) + 1;
	mut_read_seq_ss << DNATools::RandDNA(ins_size);
      }
      else if (prob - mut_prob - ins_prob <= del_prob){
	int del_size = rand()%(MAX_DEL-1) + 1;
	i += (del_size-1);
      }
      else
	mut_read_seq_ss << read_seq[i];
    }
  
    string mut_read_seq = read_seq.substr(0, PERFECT_FLANK) + mut_read_seq_ss.str() + read_seq.substr(read_seq.size()-PERFECT_FLANK, PERFECT_FLANK);

    // Align mutated read
    string ref_seq_al, read_seq_al;
    float score;
    vector<BamTools::CigarOp> cigar_list;
    NWNoRefEndPenalty::Align(ref_seq, mut_read_seq, ref_seq_al, read_seq_al, &score, cigar_list);

    if (ref_seq_al.size() != read_seq_al.size())
      PrintMessageDieOnError("Alignment lengths must match", ERROR);

    // Determine alignment start index
    unsigned int align_start = 0;
    while (align_start < read_seq_al.size() && read_seq_al[align_start] == '-')
      align_start++;
    if (align_start == read_seq_al.size()) PrintMessageDieOnError("Alignment only consists of - characters", ERROR);

    // Determine alignment end index
    int align_end = read_seq_al.size()-1;
    while (align_end >= 0 && read_seq_al[align_end] == '-')
      align_end--;
    
    int num_ins = 0;
    for (int coord = align_start; coord <= align_end; coord++)
      if (ref_seq_al[coord] == '-')
	num_ins++;
    
    // Determine if alignment boundaries are correct
    if (align_start == read_start && (align_end-num_ins-align_start+1) == READ_LEN)
      num_correct++;
    else {
      /*
      std::cerr << std::endl
		<< ref_seq_al   << std::endl
		<< read_seq_al  << std::endl 
		<< read_seq     << std::endl 
		<< mut_read_seq << std::endl
		<< align_start << " " << read_start << " " << align_end-num_ins-align_start+1 << std::endl << std::endl;
      */
    }
  }  
  return num_correct;
}



void NWNoRefEndPenaltyTest::test_Align_Mut_Only(){
  int corr_count = GenAlignments(NUM_TRIALS, 0.05, 0.0, 0.0);
  CPPUNIT_ASSERT(static_cast<float>(corr_count)/NUM_TRIALS >= MIN_FRAC_CORRECT);
}

void NWNoRefEndPenaltyTest::test_Align_Ins_Only(){
  int corr_count = GenAlignments(NUM_TRIALS, 0, 0.01, 0.0);
  CPPUNIT_ASSERT(static_cast<float>(corr_count)/NUM_TRIALS >= MIN_FRAC_CORRECT);
}

void NWNoRefEndPenaltyTest::test_Align_Del_Only(){
  int corr_count = GenAlignments(NUM_TRIALS, 0, 0, 0.01);
  CPPUNIT_ASSERT(static_cast<float>(corr_count)/NUM_TRIALS >= MIN_FRAC_CORRECT);
}

void NWNoRefEndPenaltyTest::test_Align_All_Muts(){
  int corr_count = GenAlignments(NUM_TRIALS, 0.05, 0.01, 0.01);
  CPPUNIT_ASSERT(static_cast<float>(corr_count)/NUM_TRIALS >= MIN_FRAC_CORRECT);
}
 

