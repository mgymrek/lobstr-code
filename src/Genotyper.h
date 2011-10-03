/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef GENOTYPER_H_
#define GENOTYPER_H_

#include <list>
#include <map>
#include <vector>

#include "MSReadRecord.h"
#include "MSRecord.h"
// Read with identical ReadPosition fields are assumed to be
// PCR duplicates. Note the things we are comparing all come from the 
// same MS locus, so no need to have chrom, ms start, etc.
struct ReadPosition {
  // start position of alignment
  int start;
  // length of the read
  int read_length;

  // provide overloaded operators to compare ReadPosition
  // when used as a map key
  bool operator<(const ReadPosition& read_pos1) const {
    if (start == read_pos1.start) {
      return read_length < read_pos1.read_length;
    } else {
      return start < read_pos1.start;
    }
  }
};

class Genotyper {
 public:
  Genotyper();
  ~Genotyper();

  // add an aligned read to the map
  void AddRead(MSReadRecord* read);

  // get list of copy numbers after removing pcr duplicates
  void RemovePCRDuplicates(const std::list<MSReadRecord>& str_records,
			   std::vector<float>* copy_numbers);

  // get the genotype from a list of records aligning to a 
  // single locus
  std::pair<float, float> GetGenotype(const std::vector<float>& copy_numbers,
				 bool autosomal, float ref_copy_number,
				 int* num_conflicting_reads);

  // get the genotype and locus information
  // write the genotype record to a string
  bool GetGenotypeString(const std::pair<std::string, int> & coordinate,
			 const MSReadRecord& repr_read,
			 std::string* result);

  // write genotype output to file
  void WriteOutput(const std::string& filename);

  // reset pi values
  void ResetPi(float pi0, float pi1, float pi2);

  // reset mu values
  void ResetMu(float mu0, float mu1, float mu2);

 private:
  // get probability of a genotype given read counts
  float GetPosteriorProb(const std::vector<float>& copy_numbers, int genotype,
			 float ref_copy_number);

  // keep track of aligned reads for genotyping
  // <(chrom, start), list<MSReadRecord> >
  std::map<std::pair<std::string, int>, std::list<MSReadRecord> > aligned_str_map_;

  // str loci information
  std::map<int, MSRecord> ms_loci_dict_;

  // priors
  std::map<int, float> genotype_prior;
  std::map<int, float> read_count_prior;
};

#endif /* GENOTYPER_H_ */
