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

#ifndef SRC_READCONTAINER_H_
#define SRC_READCONTAINER_H_

#include <iostream>
#include <list>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "src/cigar.h"
#include "src/ReferenceSTR.h"
#include "src/api/BamReader.h"
#include "src/api/BamMultiReader.h"

using namespace std;
using BamTools::BamReader;
using BamTools::BamMultiReader;
using BamTools::BamRegion;
using BamTools::BamAlignment;
using BamTools::SamHeader;
using BamTools::RefData;
using BamTools::RefVector;
using BamTools::CigarOp;

struct AlignedRead {
  std::string chrom;
  int msStart;
  int msEnd;
  std::string read_group; // identifies a unique sample
  int read_start;
  std::string nucleotides;
  std::string qualities;
  vector<BamTools::CigarOp> cigar_ops;
  std::string repseq;
  int period;
  int diffFromRef;
  float refCopyNum;
  int mate;
  bool strand;
  int stitched;
  int matedist;
  int mapq;
  bool stutter;
};

/*
  Class to store aligned reads from each STR locus
 */
class ReadContainer {
 public:
  ReadContainer(vector<std::string> filenames);
  ~ReadContainer();

  /* Add reads from a bam file */
  void AddReadsFromFile(const ReferenceSTR& ref_str);

  /* Get reads at an STR coordinate */
  void GetReadsAtCoord(const std::pair<std::string, int>& coord,
		       std::list<AlignedRead>* reads);

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
  int GetSTRAllele(const AlignedRead& aligned_read,
                   const CIGAR_LIST& cigar_list);

  /* Bam file reader */
  BamTools::BamMultiReader reader;
  BamTools::RefVector references;
  map<std::string, int> chrom_to_refid;
};

#endif  // SRC_READCONTAINER_H_
