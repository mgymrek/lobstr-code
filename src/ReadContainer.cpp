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

void ReadContainer::AddReadsFromFile(const ReferenceSTR& ref_str) {
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
    if (ParseRead(aln, &aligned_read)) {
      // Add to map
      pair<string, int> coord
	(aligned_read.chrom, aligned_read.msStart);
      if (aligned_str_map_.find(coord) != aligned_str_map_.end()) {
	aligned_str_map_.at(coord).push_back(aligned_read);
      } else {
	list<AlignedRead> aligned_read_list;
	aligned_read_list.push_back(aligned_read);
	aligned_str_map_.insert(pair< pair<string, int>, list<AlignedRead> >
				(coord, aligned_read_list));
      } 
    }
  }
}

bool ReadContainer::ParseRead(const BamTools::BamAlignment& aln,
			      AlignedRead* aligned_read) {
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
    if (aligned_read->mate == 0) {
      stringstream msg;
      msg << aln.Name << " from group " << aligned_read->read_group << " Could not get genotype.";
      PrintMessageDieOnError(msg.str(), ERROR);
    }
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
  cerr << aligned_read->diffFromRef << endl;
  if (aligned_read->diffFromRef + (aligned_read->refCopyNum*aligned_read->period) < MIN_ALLELE_SIZE) {
    stringstream msg;
    msg << "Discarding read " << aligned_read->ID << ". Invalid allele length";
    PrintMessageDieOnError(msg.str(), WARNING);
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
  if (aligned_str_map_.find(coord) != aligned_str_map_.end()) {
    *reads = aligned_str_map_.at(coord);
  }
}

void ReadContainer::RemovePCRDuplicates() {
  for (map<pair<string, int>, list<AlignedRead> >::iterator
         it = aligned_str_map_.begin();
       it != aligned_str_map_.end(); it++) {
    // map of <start pos, length> -> aligned reads list
    map<pair<int, bool>, list<AlignedRead> > pcr_duplicates;
    // Group into duplicates
    for (list<AlignedRead>::const_iterator
           it2 = it->second.begin(); it2 != it->second.end(); it2++) {
      //pair<int, int> key(it2->read_start, it2->nucleotides.length());
      pair<int, bool> key(it2->read_start, it2->strand);
      if (pcr_duplicates.find(key) != pcr_duplicates.end()) {
        pcr_duplicates.at(key).push_back(*it2);
      } else {
        list<AlignedRead> pcr_dup_reads;
        pcr_dup_reads.push_back(*it2);
        pcr_duplicates.insert(pair< pair<int, bool>, list<AlignedRead> >
                              (key, pcr_dup_reads));
      }
    }
    // Choose one rep from each group
    list<AlignedRead> reads_after_rmdup;
    for (map<pair<int, bool>, list<AlignedRead> >::const_iterator
           it3 = pcr_duplicates.begin();
         it3 != pcr_duplicates.end(); it3++) {
      AlignedRead rep_read;
      GetRepRead(it3->second, &rep_read);
      reads_after_rmdup.push_back(rep_read);
    }
    // Reset entry in dictionary
    aligned_str_map_.at(it->first) = reads_after_rmdup;
  }
}

void ReadContainer::GetRepRead(const list<AlignedRead>&
                               aligned_read_list,
                               AlignedRead* rep_alignment) {
  map<float, list<AlignedRead> > copy_number_to_reads;
  for (list<AlignedRead>::const_iterator
         it = aligned_read_list.begin();
       it != aligned_read_list.end(); it++) {
    float diff = it->diffFromRef;
    if (copy_number_to_reads.find(diff) !=
        copy_number_to_reads.end()) {
      copy_number_to_reads.at(diff).push_back(*it);
    } else {
      list<AlignedRead> aligned_read_list;
      aligned_read_list.push_back(*it);
      copy_number_to_reads.insert(pair<float, list<AlignedRead> >
                                  (diff, aligned_read_list));
    }
  }
  // check for majority vote, use qual as tiebreaker
  float majority_vote_copy_number = -1;
  size_t majority_vote_num_supporting_reads = 0;
  float majority_vote_average_quality = 0;
  for (map<float, list<AlignedRead> >::const_iterator
         it = copy_number_to_reads.begin();
       it != copy_number_to_reads.end(); it++) {
    float avg_qual = GetAverageQualityScore(it->second);
    if (it->second.size() > majority_vote_num_supporting_reads ||
        (it->second.size() == majority_vote_num_supporting_reads &&
         avg_qual > majority_vote_average_quality)) {
      // we are majority so far
      majority_vote_copy_number = it->first;
      majority_vote_num_supporting_reads = it->second.size();
      majority_vote_average_quality = avg_qual;
      *rep_alignment = it->second.front();
    }
  }
}

float ReadContainer::GetAverageQualityScore(const list<AlignedRead>&
                                            aligned_read_list) {
  if (aligned_read_list.size() == 0) return 0.0;
  float total_quality = 0;
  for (list<AlignedRead>::const_iterator it = aligned_read_list.begin();
       it != aligned_read_list.end(); it++) {
    total_quality += GetScore(it->qualities);
  }
  return total_quality/static_cast<float>(aligned_read_list.size());
}

float ReadContainer::GetScore(const string& quality_string) {
  if (quality_string.empty()) return 0.0;
  float total_quality = 0;
  for (size_t i = 0; i < quality_string.length(); i++) {
    total_quality +=  quality_string.at(i) - 33;
  }
  return total_quality/static_cast<float>(quality_string.length());
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
