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

#include "src/AlignedRead.h"
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

/*
  Class to store aligned reads from each STR locus
 */
class ReadContainer {
  friend class ReadContainerTest;
 public:
  ReadContainer(vector<std::string> filenames);
  ~ReadContainer();

  /* Add reads from a bam file */
  void AddReadsFromFile(const ReferenceSTR& ref_str, map<pair<string,int>, string>& ref_ext_nucleotides, const vector<string>& chroms_to_include);

  /* Clear reads from container */
  void ClearReads();

  /* Get reads at an STR coordinate */
  void GetReadsAtCoord(const std::pair<std::string, int>& coord,
		       std::list<AlignedRead>* reads);

  /* Get samples */
  void GetSampleInfo();

  // genotyper needs access to this to iterate over it
  std::map<std::pair<std::string, int>, std::list<AlignedRead> >
    aligned_str_map_;

  // list of samples
  std::vector<std::string> samples_list;
  std::map<std::string, std::string> rg_id_to_sample;

 protected:
  /* Parse BamAlignment into AlignedRead */
  bool ParseRead(const BamTools::BamAlignment& aln,
		 AlignedRead* aligned_read, 
		 map<pair<string,int>, string>& ref_ext_nucleotides);
 private:

  /* Parse bam tags into the appropriate types */
  bool GetIntBamTag(const BamTools::BamAlignment& aln,
		    const std::string& tag_name, int* destination);
  bool GetStringBamTag(const BamTools::BamAlignment& aln,
		       const std::string& tag_name, std::string* destination);
  bool GetFloatBamTag(const BamTools::BamAlignment& aln,
		      const std::string& tag_name, float* destination);

  /* Adjust diff from ref based on cigar */
  int GetSTRAllele(const CIGAR_LIST& cigar_list);

  /* Bam file reader */
  BamTools::BamMultiReader reader;
  BamTools::RefVector references;
  map<std::string, int> chrom_to_refid;
};

#endif  // SRC_READCONTAINER_H_
