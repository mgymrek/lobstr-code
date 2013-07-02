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
#include "src/VCFWriter.h"

using namespace std;

struct STRAnnotation {
  std::string chrom;
  int msStart;
  std::string name;
  std::vector<int> alleles;
};

/*
  Class to determine allelotypes at each locus
 */

class Genotyper {
 public:
  Genotyper(NoiseModel* _noise_model,
            const std::vector<std::string>& _haploid_chroms,
            std::map<pair<std::string, int>, std::string>* _ref_nucleotides,
            std::map<pair<std::string, int>, std::string>* _ref_repseq,
	    const std::string& vcf_file,
	    const std::vector<std::string>& _samples,
	    const std::map<std::string,std::string>& _rg_id_to_sample);
  ~Genotyper();

  /* Load annotations to use */
  void LoadAnnotations(const vector<std::string> annot_files);

  /* determine allelotypes and write to file */
  void Genotype(const list<AlignedRead>& read_list);

 private:
  /* Process all reads at a single locus */
  bool ProcessLocus(const std::list<AlignedRead>& aligned_reads,
                    STRRecord* str_record, bool is_haploid);

  /* Get list of alleles to use */
  bool GetAlleles(const std::list<AlignedRead>& aligned_reads,
		  std::vector<int>* alleles);

  /* Get rid of alleles with length < 0*/
  void CleanAllelesList(int reflen, std::vector<int>* alleles);

  /* Divide reads to one list per sample */
  bool GetReadsPerSample(const std::list<AlignedRead>& aligned_reads,
			 const std::vector<string>& samples,
			 const std::map<std::string,std::string>& rg_id_to_sample,
			 std::vector<std::list<AlignedRead> >* sample_reads);

  /* Get log likelihood of an allelotype */
  float CalcLogLik(int a, int b,
                const list<AlignedRead>& aligned_reads,
                int period, int* counta, int* countb);

  /* Get most likely allelotype */
  void FindMLE(const list<AlignedRead>& aligned_reads,
	       map<int,int> spanning_reads,
               bool haploid, STRRecord* str_record);

  /* chromosomes to treat as haploid */
  std::vector<std::string> haploid_chroms;

  /* store the noise model parameters */
  NoiseModel* noise_model;

  /* Reference nucleotide for each locus */
  std::map<pair<std::string, int>, std::string>* ref_nucleotides;

  /* Reference repseq for each locus */
  std::map<pair<std::string, int>, std::string>* ref_repseq;

  /* List of samples */
  std::vector<std::string> samples;
  std::map<std::string,std::string> rg_id_to_sample;

  /* annotations */
  std::map<pair<std::string, int>, STRAnnotation> annotations;

  /* File writers */
  VCFWriter* vcfWriter;
};

#endif  // SRC_GENOTYPER_H_
