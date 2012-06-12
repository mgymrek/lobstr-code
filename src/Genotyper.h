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

#ifndef SRC_GENOTYPER_H_
#define SRC_GENOTYPER_H_

#include <list>
#include <string>
#include <vector>

#include "src/NoiseModel.h"
#include "src/ReadContainer.h"

using namespace std;

/*
  Class to determine allelotypes at each locus
 */

class Genotyper {
 public:
  Genotyper(NoiseModel* _noise_model,
            bool _male, bool _simple);
  ~Genotyper();

  /* determine allelotypes and write to file */
  void Genotype(const ReadContainer& read_container,
                const std::string& output_file);

 private:
  /* Get log likelihood of an allelotype */
  float CalcLogLik(int a, int b,
                   const list<AlignedRead>& aligned_reads,
                   int period, int* counta, int* countb);

  /* Get most likely allelotype */
  void FindMLE(const list<AlignedRead>& aligned_reads, int period,
               float* allele1, float* allele2, float* score);

  /* Get allelotype without using noise model */
  void SimpleGenotype(const list<AlignedRead>& aligned_reads,
                      int period,
                      float* allele1, float* allele2, float* score);

  /* the sample is male */
  bool male;

  /* use the simple genotyper with no noise model */
  bool simple;

  /* store the noise model parameters */
  NoiseModel* noise_model;
};

#endif  // SRC_GENOTYPER_H_
