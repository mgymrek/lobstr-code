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
#include <map>
#include <set>
#include <string>
#include <vector>

#include "src/NoiseModel.h"
#include "src/ReadContainer.h"
#include "src/STRRecord.h"

using namespace std;

/*
  Class to determine allelotypes at each locus
 */

class Genotyper {
 public:
  Genotyper(NoiseModel* _noise_model,
            const std::vector<std::string>& _haploid_chroms,
            std::map<pair<std::string, int>, std::string>* _ref_nucleotides,
            std::map<pair<std::string, int>, std::string>* _ref_repseq);
  ~Genotyper();

  /* Load prior information on alleles and allele frequencies */
  void LoadPriors(const std::string& filename);

  /* determine allelotypes and write to file */
  void Genotype(const ReadContainer& read_container,
                const std::string& output_file,
                const std::string& vcf_file);

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
  void FindMLE(const list<AlignedRead>& aligned_reads,
               const map<int, float>& prior_freqs,
               bool haploid, STRRecord* str_record);

  /* Get prior for genotype <allele1,allele2>, assuming HWE */
  float GetPrior(int allele1, int allele2,
                 const map<int, float>& prior_freqs);

  /* chromosomes to treat as haploid */
  std::vector<std::string> haploid_chroms;

  /* store the noise model parameters */
  NoiseModel* noise_model;

  /* Use information about known alleles */
  bool use_known_alleles;

  /* Information on allele frequencies to use as priors */
  map<pair<string, int>, map<int, float> > allele_frequencies_per_locus;

  /* Reference nucleotide for each locus */
  std::map<pair<std::string, int>, std::string>* ref_nucleotides;

  /* Reference repseq for each locus */
  std::map<pair<std::string, int>, std::string>* ref_repseq;
};

#endif  // SRC_GENOTYPER_H_
