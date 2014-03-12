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

#include <err.h>

#include <iostream>
#include <list>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "src/common.h"
#include "src/ReadContainer.h"
#include "src/AlignmentFilters.h"
#include "src/runtime_parameters.h"

using namespace std;

const int MIN_ALLELE_SIZE = 0;

ReadContainer::ReadContainer(vector<std::string> filenames) {
  string bamfile;
  vector<string> index_files;
  // Open bam files
  for (size_t i = 0; i < filenames.size(); i++) {
    bamfile = filenames.at(i);
    index_files.push_back(bamfile + ".bai");
    if (!reader.OpenFile(bamfile)) {
      PrintMessageDieOnError("Could not open " + bamfile, ERROR);
    }
  }
  // Check indexes
  if (!reader.OpenIndexes(index_files)) {
    PrintMessageDieOnError("Could not open index files", ERROR);
  }
  // get chrom_to_refid
  references = reader.GetReferenceData();
  for (size_t i = 0; i < references.size(); i++) {
    chrom_to_refid[references.at(i).RefName] = static_cast<int>(i);
  }
}

void ReadContainer::AddReadsFromFile(const ReferenceSTR& ref_str, map<pair<string,int>, string>& ref_ext_nucleotides) {
  if (ref_str.chrom != "NA") {
    int refid = -1;
    if (chrom_to_refid.find(ref_str.chrom) !=
	chrom_to_refid.end()) {
      refid = chrom_to_refid.at(ref_str.chrom);
    }
    if (refid == -1) {
      PrintMessageDieOnError("Could not locate STR reference chromosome in bam file", ERROR);
    }
    BamTools::BamRegion bam_region(refid, ref_str.start-extend, refid, ref_str.stop+extend);
    if (!reader.SetRegion(bam_region)) {
      PrintMessageDieOnError("Could not set bam region", ERROR);
    }
  }
  BamTools::BamAlignment aln;
  while (reader.GetNextAlignment(aln)) {
    AlignedRead aligned_read;
    if (ParseRead(aln, &aligned_read, ref_ext_nucleotides)) {
      // Add to map
      pair<string, int> coord
	(aligned_read.chrom, aligned_read.msStart);
      if (aligned_str_map_.find(coord) != aligned_str_map_.end()) {
	aligned_str_map_.at(coord).push_back(aligned_read);
      } 
      else {
	list<AlignedRead> aligned_read_list;
	aligned_read_list.push_back(aligned_read);
	aligned_str_map_.insert(pair< pair<string, int>, list<AlignedRead> >
				(coord, aligned_read_list));
      } 
    }
  }
}

bool ReadContainer::ParseRead(const BamTools::BamAlignment& aln,
			      AlignedRead* aligned_read, 
			      map<pair<string,int>, string>& ref_ext_nucleotides) {
  // get read ID
  aligned_read->ID = aln.Name;
  // get nucleotides
  aligned_read->nucleotides = aln.QueryBases;
  // get qualities
  aligned_read->qualities = aln.Qualities;
  // get strand
  aligned_read->strand = aln.IsReverseStrand();
  // get chrom
  aligned_read->chrom = references.at(aln.RefID).RefName;
  // get read start
  aligned_read->read_start = aln.Position;
  // get cigar
  aligned_read->cigar_ops = aln.CigarData;
  // get if mate pair
  if (aln.IsSecondMate()) {
    aligned_read->mate = 1;
  } else {
    aligned_read->mate = 0;
  }
  // Only process if it is the primary alignment
  if (aligned_read->mate) {
    return false;
  }
  // Get all the tag data
  // don't process if partially spanning (from old lobSTR)
  int partial = 0;
  if (GetIntBamTag(aln, "XP", &partial)) {
    if (partial == 1) return false;
  }
  // get read group
  if (!GetStringBamTag(aln, "RG", &aligned_read->read_group)) {
    stringstream msg;
    msg << aln.Name << " Could not get read group.";
    PrintMessageDieOnError(msg.str(), ERROR);
  }
  // get msStart
  if (!GetIntBamTag(aln, "XS", &aligned_read->msStart)) {
    stringstream msg;
    msg << aln.Name << " from group " << aligned_read->read_group << " Could not get STR start coordinate. Did this bam file come from lobSTR?";
    PrintMessageDieOnError(msg.str(), ERROR);
  }
  // get msEnd
  if (!GetIntBamTag(aln, "XE", &aligned_read->msEnd)) {
    stringstream msg;
    msg << aln.Name << " from group " << aligned_read->read_group << " Could not get STR end coordinate. Did this bam file come from lobSTR?";
    PrintMessageDieOnError(msg.str(), ERROR);
  }
  // get mapq. Try unsigned/signed
  if (!GetIntBamTag(aln, "XQ", &aligned_read->mapq)) {
    stringstream msg;
    aligned_read->mapq = 0;
  }
  // get diff
  if (!GetIntBamTag(aln, "XD", &aligned_read->diffFromRef)) {
    return false;
  }
  // get mate dist
  if (!GetIntBamTag(aln, "XM", &aligned_read->matedist)) {
    aligned_read->matedist = 0;
  }
  // get STR seq
  if (!GetStringBamTag(aln, "XR", &aligned_read->repseq)) {
    stringstream msg;
    msg << aln.Name << " from group " << aligned_read->read_group << " Could not get repseq.";
    PrintMessageDieOnError(msg.str(), ERROR);
  }
  // get if stitched
  if (!GetIntBamTag(aln, "XX", &aligned_read->stitched)) {
    aligned_read->stitched = 0;
  }
  // get ref copy num
  if (!GetFloatBamTag(aln, "XC", &aligned_read->refCopyNum)) {
    stringstream msg;
    msg << aln.Name << " from group " << aligned_read->read_group << " Could not get reference copy number.";
    PrintMessageDieOnError(msg.str(), ERROR);
  }
  // get period
  aligned_read->period = aligned_read->repseq.length();
  if (include_flank) {  // diff is just sum of differences in cigar
    CIGAR_LIST cigar_list;
    for (vector<BamTools::CigarOp>::const_iterator
	   it = aligned_read->cigar_ops.begin();
	 it != aligned_read->cigar_ops.end(); it++) {
      CIGAR cig;
      cig.num = (*it).Length;
      cig.cigar_type = (*it).Type;
      cigar_list.cigars.push_back(cig);
    }
    bool added_s;
    bool cigar_had_s;
    cigar_list.ResetString();
    GenerateCorrectCigar(&cigar_list, aln.QueryBases,
			 &added_s, &cigar_had_s);
    aligned_read->diffFromRef = GetSTRAllele(cigar_list);
  }
  // apply filters
  if (unit) {
    if (aligned_read->diffFromRef % aligned_read->period != 0) return false;
  }
  if (abs(aligned_read->diffFromRef) > max_diff_ref) {
    return false;
  }
  if (aligned_read->mapq > max_mapq) {
    return false;
  }
  if (aligned_read->matedist > max_matedist) {
    return false;
  }
  // Check if the allele length is valid
  if (aligned_read->diffFromRef + (aligned_read->refCopyNum*aligned_read->period) < MIN_ALLELE_SIZE) {
    return false;
  }

  // check that read sufficiently spans STR
  int max_read_start = aligned_read->msStart - min_border;
  int min_read_stop  = aligned_read->msEnd   + min_border;
  if (aln.Position > max_read_start || aln.GetEndPosition() < min_read_stop)
    return false;
  
  // check that both ends of the read contain sufficient perfect matches
  if (min_read_end_match > 0){
    map<pair<string,int>, string>::iterator loc_iter = ref_ext_nucleotides.find(pair<string,int>(aligned_read->chrom, aligned_read->msStart));
    if (loc_iter == ref_ext_nucleotides.end())
      PrintMessageDieOnError("No extended reference sequence found for locus", ERROR);
    
    string ref_ext_seq = loc_iter->second;
    pair<int,int> num_end_matches = GetNumEndMatches(aligned_read, ref_ext_seq, aligned_read->msStart-extend);
    if (num_end_matches.first < min_read_end_match || num_end_matches.second < min_read_end_match)
      return false;
  }
  return true;
}

bool ReadContainer::GetIntBamTag(const BamTools::BamAlignment& aln,
		  const std::string& tag_name, int* destination) {
  char tag_type;
  if (!aln.GetTagType(tag_name, tag_type)) {return false;}
  switch (tag_type) {
  case (BamTools::Constants::BAM_TAG_TYPE_INT32):
    return aln.GetTag(tag_name, *destination);
  case (BamTools::Constants::BAM_TAG_TYPE_INT8):
    int8_t d8;
    if (!aln.GetTag(tag_name, d8)) {
      return false;
    }
    *destination = static_cast<int>(d8);
    return true;
  case (BamTools::Constants::BAM_TAG_TYPE_UINT8):
    uint8_t ud8;
    if (!aln.GetTag(tag_name, ud8)) {
      return false;
    }
    *destination = static_cast<int>(ud8);
    return true;
  case (BamTools::Constants::BAM_TAG_TYPE_INT16):
    int16_t d16;
    if (!aln.GetTag(tag_name, d16)) {
      return false;
    }
    *destination = static_cast<int>(d16);
    return true;
  case (BamTools::Constants::BAM_TAG_TYPE_UINT16):
    uint16_t ud16;
    if (!aln.GetTag(tag_name, ud16)) {
      return false;
    }
    *destination = static_cast<int>(ud16);
    return true;
  case (BamTools::Constants::BAM_TAG_TYPE_UINT32):
    uint32_t ud32;
    if (!aln.GetTag(tag_name, ud32)) {
      return false;
    }
    *destination = static_cast<int>(ud32);
    return true;
  default:
    stringstream msg;
    msg << "Encountered unsupported tag type " << tag_type;
    PrintMessageDieOnError(msg.str(), ERROR);
  }
  return false;
}

bool ReadContainer::GetStringBamTag(const BamTools::BamAlignment& aln,
		     const std::string& tag_name, std::string* destination) {
  if (!aln.GetTag(tag_name, *destination)) {
    return false;
  }
  return true;
}

bool ReadContainer::GetFloatBamTag(const BamTools::BamAlignment& aln,
		     const std::string& tag_name, float* destination) {
  if (!aln.GetTag(tag_name, *destination)) {
    return false;
  }
  return true;
}

void ReadContainer::ClearReads() {
  aligned_str_map_.clear();
}

void ReadContainer::GetReadsAtCoord(const pair<string,int>& coord,
				    list<AlignedRead>* reads) {
  reads->clear();
  if (aligned_str_map_.find(coord) != aligned_str_map_.end()) {
    *reads = aligned_str_map_.at(coord);
  }
}

int ReadContainer::GetSTRAllele(const CIGAR_LIST& cigar_list) {
  int diff_from_ref = 0;
  for (size_t i = 0; i < cigar_list.cigars.size(); i++) {
    if (cigar_list.cigars.at(i).cigar_type == 'I') {
      diff_from_ref += cigar_list.cigars.at(i).num;
    }
    if (cigar_list.cigars.at(i).cigar_type == 'D') {
      diff_from_ref -= cigar_list.cigars.at(i).num;
    }
  }
  // set STR region
  return diff_from_ref;
}

ReadContainer::~ReadContainer() {}
