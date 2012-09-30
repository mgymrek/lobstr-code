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

#ifndef SRC_GENOTYPER_V2_H_
#define SRC_GENOTYPER_V2_H_

#include <list>
#include <string>
#include <vector>

#include "src/NoiseModel.h"
#include "src/ReadContainer.h"

using namespace std;

/*
  Struct to keep track of info for a single STR locus
 */
struct STRRecord {
  std::string chrom;
  int start;
  int stop;
  std::string repseq;
  int period;
  float allele1;
  float allele2;
  int coverage;
  float score;
  float allele1_score;
  float allele2_score;
  int conflicting;
  int agreeing;
  int partial_coverage;
  int num_stitched;
  float refcopy;
  int max_partial;
  std::string readstring;
  std::string partialreadstring;
  std::string max_partial_string;
  std::string allele1_string;
  std::string allele2_string;
  void Reset() {
    chrom = "";
    start = -1;
    stop = -1;
    repseq = "";
    period = -1;
    allele1 = -10000;
    allele2 = -10000;
    coverage = 0;
    allele1_score = -1;
    allele2_score = -1;
    conflicting = 0;
    agreeing = 0;
    partial_coverage = 0;
    num_stitched = 0;
    refcopy = -1;
    max_partial = 0;
  }
};

/*
  Class to determine allelotypes at each locus
 */

class Genotyper {
 public:
  Genotyper(NoiseModel* _noise_model,
            const std::vector<std::string>& _haploid_chroms,
            bool _simple);
  ~Genotyper();

  /* determine allelotypes and write to file */
  void Genotype(const ReadContainer& read_container,
                const std::string& output_file);

 private:
  /* Process all reads at a single locus */
  bool ProcessLocus(const std::string chrom,
                    const int str_coord,
                    const std::list<AlignedRead>& aligned_reads,
                    STRRecord* str_record);

  /* Get log likelihood of an allelotype */
  float CalcLogLik(int a, int b,
                const list<AlignedRead>& aligned_reads,
                int period, int* counta, int* countb);

  /* Get most likely allelotype */
  void FindMLE(const list<AlignedRead>& aligned_reads, int period,
               float* allele1, float* allele2, float* score,
               float* score_allele1, float* score_allele2, bool haploid);

  /* Get allelotype without using noise model */
  void SimpleGenotype(const list<AlignedRead>& aligned_reads,
                      int period,
                      float* allele1, float* allele2, float* score);

  /* chromosomes to treat as haploid */
  std::vector<std::string> haploid_chroms;

  /* use the simple genotyper with no noise model */
  bool simple;

  /* store the noise model parameters */
  NoiseModel* noise_model;
};

#endif  // SRC_GENOTYPER_V2_H_
