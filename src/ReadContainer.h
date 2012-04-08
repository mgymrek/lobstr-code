/*
 * Author: Melissa Gymrek 2012
 */

#ifndef READ_CONTAINER_H_
#define READ_CONTAINER_H_

#include <iostream>
#include <list>
#include "cigar.h"
#include "api/BamReader.h"

using namespace std;
using BamTools::BamReader;
using BamTools::BamAlignment;
using BamTools::SamHeader;
using BamTools::RefData;
using BamTools::RefVector;
using BamTools::CigarOp;
struct AlignedRead {
  std::string chrom;
  int msStart;
  int msEnd;
  int read_start;
  std::string nucleotides;
  std::string qualities;
  vector<BamTools::CigarOp> cigar_ops;
  std::string repseq;
  int period;
  int diffFromRef;
  float refCopyNum;
  int partial;
  bool strand;
};

/*
  Class to store aligned reads from each STR locus
 */
class ReadContainer {
 public:
  ReadContainer();
  ~ReadContainer();

  /* Add reads from a bam file */
  void AddReadsFromFile(std::string bamfile);

  /* Remove pcr duplicates */
  void RemovePCRDuplicates();

  // genotyper needs access to this to iterate over it
  std::map<std::pair<std::string, int>, std::list<AlignedRead> > 
    aligned_str_map_;

 private:
  /* Get values from representative read in set of dups */
  void GetRepRead(const list<AlignedRead>& aligned_read_list,
		  AlignedRead* rep_alignment);

  /* Get average quality score of a set of reads */
  float GetAverageQualityScore(const list<AlignedRead>&
			       aligned_read_list);

  /* Get quality core for a single read */
  float GetScore(const std::string& quality_string);

  /* Adjust diff from ref based on cigar */
  int GetSTRAllele(const AlignedRead& aligned_read, const CIGAR_LIST& cigar_list);

  BamTools::BamReader reader;
};

#endif /* READ_CONTAINER_H_ */
