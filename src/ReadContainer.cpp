/*
Copyright (C) 2011-2014 Melissa Gymrek <mgymrek@mit.edu>
Revisions     2014 Thomas Willems <twillems@mit.edu> 

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

#include <algorithm>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "src/common.h"
#include "src/nw.h"
#include "src/ReadContainer.h"
#include "src/AlignmentFilters.h"
#include "src/runtime_parameters.h"

using namespace std;

const int MIN_ALLELE_SIZE = 0;
const int CIGAR_BUFFER = 5;

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
  // Get sample info
  GetSampleInfo();
  // Open writers
  if (output_bams) {
    map<string, int> chrom_sizes;
    for (size_t i = 0; i < references.size(); i++) {
      chrom_sizes[references.at(i).RefName] = (int)references.at(i).RefLength;
    }
    writer_reads = new SamFileWriter(output_prefix + ".reads.bam", chrom_sizes);
    writer_filtered = new SamFileWriter(output_prefix + ".filtered.bam", chrom_sizes);
  }
}

void ReadContainer::GetSampleInfo() {
  BamTools::SamHeader header = reader.GetHeader();
  if (!header.HasReadGroups()) {
    PrintMessageDieOnError("No read groups found", ERROR);
  }
  BamTools::SamReadGroupDictionary read_groups = header.ReadGroups;
  for (BamTools::SamReadGroupIterator it = read_groups.Begin();
       it != read_groups.End(); it++) {
    const BamTools::SamReadGroup& rg = *it;
    string rg_sample;
    if (!rg.HasSample()) {
      PrintMessageDieOnError("No sample in read group for " + rg.ID, WARNING);
      rg_sample = rg.ID;
    } else {
      rg_sample = rg.Sample;
    }
    if (my_verbose) {
      PrintMessageDieOnError("Adding sample " + rg_sample + " " + rg.ID, PROGRESS);
    }
    rg_id_to_sample.insert(pair<string,string>(rg.ID,rg_sample));
    if (find(samples_list.begin(), samples_list.end(), rg_sample) ==
    	samples_list.end()) {
      samples_list.push_back(rg_sample);
    }
  }
}

void ReadContainer::AddReadsFromFile(const ReferenceSTR& ref_str, const vector<ReferenceSTR>& ref_str_chunk,
                                     map<pair<string,int>, string>& ref_ext_nucleotides,
                                     const vector<string>& chroms_to_include) {
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
  STRIntervalTree itree;
  itree.LoadIntervals(ref_str_chunk);
  while (reader.GetNextAlignment(aln)) {
    if (chroms_to_include.size() > 0 && find(chroms_to_include.begin(), chroms_to_include.end(), references.at(aln.RefID).RefName) == chroms_to_include.end()) {
      continue;
    }
    vector<AlignedRead> aligned_reads;
    if (!ParseRead(aln, &aligned_reads, itree, ref_ext_nucleotides)) {
      continue;
    }
    for (vector<AlignedRead>::const_iterator it = aligned_reads.begin();
         it != aligned_reads.end(); it++) {
      const AlignedRead& aligned_read = *it;
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
                              vector<AlignedRead>* aligned_reads,
                              STRIntervalTree& itree,
                              map<pair<string,int>, string>& ref_ext_nucleotides) {
  // Dummy aligned read to set fields common to all
  AlignedRead dummy_aligned_read;
  // Keep track of STRs already covered so don't add twice (for annotated markers with two listings
  set<int> str_starts;
  // *** First check bam flags *** //
  if (!aln.IsMapped()) {
    if (output_bams) {
      writer_filtered->WriteAllelotypeRead(aln, "UNMAPPED", "", -1, -1, "", 0);
    }
    return false;
  }
  // *** Get read properties*** //
  // get read ID
  dummy_aligned_read.ID = aln.Name;
  // get nucleotides
  dummy_aligned_read.nucleotides = aln.QueryBases;
  // get qualities
  dummy_aligned_read.qualities = aln.Qualities;
  // get strand
  dummy_aligned_read.strand = aln.IsReverseStrand();
  // get chrom
  dummy_aligned_read.chrom = references.at(aln.RefID).RefName;
  // get read start
  dummy_aligned_read.read_start = aln.Position;
  // get cigar
  dummy_aligned_read.cigar_ops = aln.CigarData;
  CIGAR_LIST cigar_list;
  if (!GetCigarList(dummy_aligned_read, &cigar_list)) {
    if (output_bams) {
      writer_filtered->WriteAllelotypeRead(aln, "NO_CIGAR", dummy_aligned_read.chrom, -1, -1, "", 0);
    }
    return false;
  }
  // Ignore pair information
  dummy_aligned_read.mate = 0;
  // get mapq
  if (!GetIntBamTag(aln, "XQ", &dummy_aligned_read.mapq)) {
    dummy_aligned_read.mapq = 0;
  }
  // get mate dist
  if (!GetIntBamTag(aln, "XM", &dummy_aligned_read.matedist)) {
    dummy_aligned_read.matedist = 0;
  }
  // Ignore stitch information
  if (!GetIntBamTag(aln, "XX", &dummy_aligned_read.stitched)) {
    dummy_aligned_read.stitched = 0;
  }
  // get read group
  if (!GetStringBamTag(aln, "RG", &dummy_aligned_read.read_group)) {
    stringstream msg;
    msg << aln.Name << " Could not get read group.";
    PrintMessageDieOnError(msg.str(), ERROR);
  }
  // *** Alignment filters (these don't depend on which STR aligned to) *** //
  if (dummy_aligned_read.mapq > max_mapq) {
    filter_counter.increment(FilterCounter::MAPPING_QUALITY);
    if (output_bams) {
      writer_filtered->WriteAllelotypeRead(aln, filter_counter.GetFilterType(FilterCounter::MAPPING_QUALITY), dummy_aligned_read.chrom, -1, -1, "", 0);
    }
    return false;
  }
  if (aln.MapQuality == 0 && filter_mapq0) {
    filter_counter.increment(FilterCounter::MAPQ0);
    if (output_bams) {
      writer_filtered->WriteAllelotypeRead(aln, filter_counter.GetFilterType(FilterCounter::MAPQ0), dummy_aligned_read.chrom, -1, -1, "", 0);
    }
    return false;
  }
  if (dummy_aligned_read.matedist > max_matedist) {
    filter_counter.increment(FilterCounter::MATE_DIST);
    if (output_bams) {
      writer_filtered->WriteAllelotypeRead(aln, filter_counter.GetFilterType(FilterCounter::MATE_DIST), dummy_aligned_read.chrom, -1, -1, "", 0);
    }
    return false;
  }
  if (dummy_aligned_read.nucleotides.find("N") != std::string::npos) { 
    if (filter_reads_with_n) {
      filter_counter.increment(FilterCounter::CONTAINS_N_BASE);
      if (output_bams) {
        writer_filtered->WriteAllelotypeRead(aln, filter_counter.GetFilterType(FilterCounter::CONTAINS_N_BASE), dummy_aligned_read.chrom, -1, -1, "", 0);
      }
      return false;
    }
  }
  if (cigar_list.cigar_string.find("S") != std::string::npos ||
      cigar_list.cigar_string.find("H") != std::string::npos) {
    if (filter_clipped) {
      filter_counter.increment(FilterCounter::CONTAINS_CLIP);
      if (output_bams) {
        writer_filtered->WriteAllelotypeRead(aln, filter_counter.GetFilterType(FilterCounter::CONTAINS_CLIP), dummy_aligned_read.chrom, -1, -1, "", 0);
      }
      return false;
    }
  }
  // check that both ends of the aligned read have sufficient bases before the first indel
  if (min_bp_before_indel > 0){
    pair<int, int> num_bps = AlignmentFilters::GetEndDistToIndel(&dummy_aligned_read);
    if (num_bps.first != -1 && num_bps.first < min_bp_before_indel){
      filter_counter.increment(FilterCounter::BP_BEFORE_INDEL);
      if (output_bams) {
        writer_filtered->WriteAllelotypeRead(aln, filter_counter.GetFilterType(FilterCounter::BP_BEFORE_INDEL), dummy_aligned_read.chrom, -1, -1, "", 0);
      }
      return false;
    }
    if (num_bps.second != -1 && num_bps.second < min_bp_before_indel){
      filter_counter.increment(FilterCounter::BP_BEFORE_INDEL);
      if (output_bams) {
        writer_filtered->WriteAllelotypeRead(aln, filter_counter.GetFilterType(FilterCounter::BP_BEFORE_INDEL), dummy_aligned_read.chrom, -1, -1, "", 0);
      }
      return false;
    }
  }
  // *** Determine region spanned by this read *** //
  int read_start = dummy_aligned_read.read_start;
  int read_end = dummy_aligned_read.read_start + (int)(dummy_aligned_read.nucleotides.size()) - GetSTRAllele(cigar_list);
  // *** Determine which reference STRs overlapped by this read *** //
  vector<ReferenceSTR> spanned_strs;
  itree.GetSpannedIntervals(read_start, read_end, &spanned_strs);
  if (spanned_strs.size() == 0)  {
    if (output_bams) {
      writer_filtered->WriteAllelotypeRead(aln, "NO_SPANNED_STRS", dummy_aligned_read.chrom, -1, -1, "", 0);
    }
    return false;
  }
  // realign before proceeding if specified. Do here because we need the start coord
  if (realign) {
    dummy_aligned_read.msStart = spanned_strs.at(0).start;
    if (!RedoLocalAlignment(&dummy_aligned_read, ref_ext_nucleotides)) {
      stringstream msg;
      msg << "Realignment of " << dummy_aligned_read.ID << " failed";
      PrintMessageDieOnError(msg.str(), WARNING);
      return false;
    }
    cigar_list.Clear();
    if (!GetCigarList(dummy_aligned_read, &cigar_list)) {
      if (output_bams) {
        writer_filtered->WriteAllelotypeRead(aln, "NO_CIGAR", dummy_aligned_read.chrom, -1, -1, "", 0);
      }
      return false;
    }
  }
  for (size_t i = 0; i < spanned_strs.size(); i++) {
    const ReferenceSTR ref_str = spanned_strs.at(i);
    AlignedRead aligned_read = dummy_aligned_read; // copy what's in dummy_aligned_read
    // get msStart
    aligned_read.msStart = ref_str.start;
    // get msEnd
    aligned_read.msEnd = ref_str.stop;
    // get STR seq
    aligned_read.repseq = ref_str.motif;
    // get ref copy num
    aligned_read.refCopyNum = static_cast<float>(ref_str.stop - ref_str.start + 1)/static_cast<float>(ref_str.motif.size());
    // get allele
    CIGAR_LIST str_cigar_list;
    if (!ExtractCigar(cigar_list, aligned_read.read_start+1, ref_str.start-CIGAR_BUFFER, ref_str.stop+CIGAR_BUFFER, &str_cigar_list)) {
      if (output_bams) {
        writer_filtered->WriteAllelotypeRead(aln, "ERROR_EXTRACTING_CIGAR", ref_str.chrom, ref_str.start, ref_str.stop, ref_str.motif, 0);
      }
      continue;
    }
    aligned_read.diffFromRef = GetSTRAllele(str_cigar_list);
    // get period
    aligned_read.period = aligned_read.repseq.length(); 
    // apply filters
    if (unit) {
      if (aligned_read.diffFromRef % aligned_read.period != 0){ 
        filter_counter.increment(FilterCounter::NOT_UNIT);
        if (output_bams) {
          writer_filtered->WriteAllelotypeRead(aln, filter_counter.GetFilterType(FilterCounter::NOT_UNIT),
                                               ref_str.chrom, ref_str.start, ref_str.stop, ref_str.motif, 0);
        }
        continue;
      }
    }
    if (abs(aligned_read.diffFromRef) > max_diff_ref) {
      filter_counter.increment(FilterCounter::DIFF_FROM_REF);
      if (output_bams) {
        writer_filtered->WriteAllelotypeRead(aln, filter_counter.GetFilterType(FilterCounter::DIFF_FROM_REF),
                                             ref_str.chrom, ref_str.start, ref_str.stop, ref_str.motif, 0);
      }
      continue;
    }
    // Check if the allele length is valid
    if (aligned_read.diffFromRef + (aligned_read.refCopyNum*aligned_read.period) < MIN_ALLELE_SIZE) {
      filter_counter.increment(FilterCounter::ALLELE_SIZE);
      if (output_bams) {
        writer_filtered->WriteAllelotypeRead(aln, filter_counter.GetFilterType(FilterCounter::ALLELE_SIZE),
                                             ref_str.chrom, ref_str.start, ref_str.stop, ref_str.motif, 0);
      }
      continue;
    }
    // check that read sufficiently spans STR
    int max_read_start = aligned_read.msStart - min_border;
    int min_read_stop  = aligned_read.msEnd   + min_border;
    if (aln.Position > max_read_start || aln.GetEndPosition() < min_read_stop){
      filter_counter.increment(FilterCounter::SPANNING_AMOUNT);
      if (output_bams) {
        writer_filtered->WriteAllelotypeRead(aln, filter_counter.GetFilterType(FilterCounter::SPANNING_AMOUNT),
                                             ref_str.chrom, ref_str.start, ref_str.stop, ref_str.motif, 0);
      }
      continue;
    }
    // Check that the read ends don't contain too many repeat instances
    if (max_repeats_in_ends > -1) {
      if (AlignmentFilters::GetMaxRepeatsInEnds(&aligned_read, aligned_read.repseq.length()*4) > 
          max_repeats_in_ends) {
        filter_counter.increment(FilterCounter::MAX_REPEATS_IN_ENDS);
        if (output_bams) {
          writer_filtered->WriteAllelotypeRead(aln, filter_counter.GetFilterType(FilterCounter::MAX_REPEATS_IN_ENDS),
                                               ref_str.chrom, ref_str.start, ref_str.stop, ref_str.motif, 0);
        }
        continue;
      }
    }
    if (str_starts.find(aligned_read.msStart) == str_starts.end()) {
      aligned_reads->push_back(aligned_read);
      str_starts.insert(aligned_read.msStart);
    }
  }
  if (aligned_reads->size() == 0) {
    if (output_bams) {
      writer_filtered->WriteAllelotypeRead(aln, "NO_SPANNED_STRS", dummy_aligned_read.chrom, -1, -1, "", 0);
    }
    return false;
  }
  // *** Apply rest of read-level filters. Get sequence based on first aligned read location *** //
  // check that the prefix and suffix of the read match maximally compared to proximal reference locations
  if (maximal_end_match_window > 0){
    map<pair<string,int>, string>::iterator loc_iter = 
      ref_ext_nucleotides.find(pair<string,int>(aligned_reads->at(0).chrom,
						aligned_reads->at(0).msStart));
    if (loc_iter == ref_ext_nucleotides.end()) {
      stringstream msg;
      msg << "No extended reference sequence found for locus " << aligned_reads->at(0).chrom << ":"
          << aligned_reads->at(0).msStart << " read " << dummy_aligned_read.ID;
      PrintMessageDieOnError(msg.str(), WARNING);
      if (output_bams) {
        writer_filtered->WriteAllelotypeRead(aln, "NO_REF_SEQ", dummy_aligned_read.chrom, -1, -1, "", 0);
      }
      return false;
    }
    string ref_ext_seq = loc_iter->second;
    bool maximum_end_matches = AlignmentFilters::HasLargestEndMatches(&dummy_aligned_read, ref_ext_seq, aligned_reads->at(0).msStart-extend,
								      maximal_end_match_window, maximal_end_match_window);
    if (!maximum_end_matches){
      filter_counter.increment(FilterCounter::NOT_MAXIMAL_END);
      if (output_bams) {
        writer_filtered->WriteAllelotypeRead(aln, filter_counter.GetFilterType(FilterCounter::NOT_MAXIMAL_END), dummy_aligned_read.chrom, -1, -1, "", 0);
      }
      return false;
    }
  }
  // check that both ends of the read contain sufficient perfect matches
  if (min_read_end_match > 0){
    map<pair<string,int>, string>::iterator loc_iter = 
      ref_ext_nucleotides.find(pair<string,int>(aligned_reads->at(0).chrom, aligned_reads->at(0).msStart));
    if (loc_iter == ref_ext_nucleotides.end()) {
      stringstream msg;
      msg << "No extended reference sequence found for locus " << aligned_reads->at(0).chrom << ":"
          << aligned_reads->at(0).msStart << " read " << aligned_reads->at(0).ID;
      PrintMessageDieOnError(msg.str(), WARNING);
      if (output_bams) {
        writer_filtered->WriteAllelotypeRead(aln, "NO_REF_SEQ", dummy_aligned_read.chrom, -1, -1, "", 0);
      }
      return false;
    }
    string ref_ext_seq = loc_iter->second;
    pair<int,int> num_end_matches = AlignmentFilters::GetNumEndMatches(&dummy_aligned_read, ref_ext_seq, aligned_reads->at(0).msStart-extend);
    if (num_end_matches.first < min_read_end_match || num_end_matches.second < min_read_end_match){
      filter_counter.increment(FilterCounter::NUM_END_MATCHES);
      if (output_bams) {
        writer_filtered->WriteAllelotypeRead(aln, filter_counter.GetFilterType(FilterCounter::NUM_END_MATCHES), dummy_aligned_read.chrom, -1, -1, "", 0);
      }
      return false;
    }
  }
  
  if (aligned_reads->size() > 0) {
    filter_counter.increment(FilterCounter::UNFILTERED);
    for (size_t i = 0; i < aligned_reads->size(); i++) {
      if (output_bams) {
        writer_reads->WriteAllelotypeRead(aln, "PASS", aligned_reads->at(i).chrom,
                                          aligned_reads->at(i).msStart, aligned_reads->at(i).msEnd,
                                          aligned_reads->at(i).repseq, aligned_reads->at(i).diffFromRef);
      }
    }
    return true;
  }
  return false;
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

bool ReadContainer::GetCigarList(const AlignedRead& aligned_read,
				 CIGAR_LIST* cigar_list) {
  for (vector<BamTools::CigarOp>::const_iterator
	 it = aligned_read.cigar_ops.begin();
       it != aligned_read.cigar_ops.end(); it++) {
    CIGAR cig;
    cig.num = (*it).Length;
    cig.cigar_type = (*it).Type;
    cigar_list->cigars.push_back(cig);
  }
  bool added_s, cigar_had_s;
  cigar_list->ResetString();
  GenerateCorrectCigar(cigar_list, aligned_read.nucleotides, &added_s, &cigar_had_s);
  return true;
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

bool ReadContainer::RedoLocalAlignment(AlignedRead* aligned_read,
                                       const map<pair<string, int>, string>& ref_ext_nucleotides) {
  // Check that we have reference sequence
  if (ref_ext_nucleotides.find(pair<string, int>(aligned_read->chrom, aligned_read->msStart)) ==
      ref_ext_nucleotides.end()) {
    stringstream msg;
    msg << "No extended reference sequence found for locus "<< aligned_read->chrom << ":"
        << aligned_read->msStart << " read " << aligned_read->ID;
    PrintMessageDieOnError(msg.str(), WARNING);
    return false;
  }
  // Pull out reference sequence
  const std::string& ref_ext_seq = ref_ext_nucleotides.at(pair<string,int>(aligned_read->chrom,
                                                                           aligned_read->msStart));
  CIGAR_LIST read_cigar_list, cigar_list;
  if (!GetCigarList(*aligned_read, &read_cigar_list)) {
    return false;
  }
  const int pad=50;
  const int read_span = int(aligned_read->nucleotides.size()) - GetSTRAllele(read_cigar_list);
  int ref_seq_start = aligned_read->msStart - extend;
  int ref_index = aligned_read->read_start - ref_seq_start;
  const std::string ref_seq = ref_ext_seq.substr(ref_index-pad, read_span+2*pad);
  // Get aligned sequence
  const std::string aligned_seq = aligned_read->nucleotides;
  // Run local realignment
  string aligned_seq_sw, ref_seq_sw;
  int sw_score;
  nw(aligned_seq, ref_seq, aligned_seq_sw, ref_seq_sw,
     &sw_score, &cigar_list);
  cigar_list.ResetString();
  // get rid of end gaps and update coords
  aligned_read->read_start = aligned_read->read_start - pad;
  if (cigar_list.cigars.at(0).cigar_type == 'D') {
    const int& num = cigar_list.cigars.at(0).num;
    aligned_read->read_start += num;
    cigar_list.cigars.erase(cigar_list.cigars.begin());
  }
  if (cigar_list.cigars.at(cigar_list.cigars.size() - 1).cigar_type == 'D') {
    cigar_list.cigars.erase(cigar_list.cigars.end() - 1);
  }
  if (cigar_list.cigars.at(0).cigar_type == 'I') {
    cigar_list.cigars.at(0).cigar_type = 'S';
  }
  if (cigar_list.cigars.at(cigar_list.cigars.size() - 1).cigar_type == 'D') {
    cigar_list.cigars.at(cigar_list.cigars.size() - 1).cigar_type = 'S';
  }
  cigar_list.ResetString();
  // make sure CIGAR is valid
  bool added_s;
  bool cigar_had_s;
  GenerateCorrectCigar(&cigar_list, aligned_read->nucleotides, &added_s, &cigar_had_s);
  vector<BamTools::CigarOp> cigar_data;
  for (size_t i = 0; i < cigar_list.cigars.size(); i++) {
    char cigar_type = cigar_list.cigars.at(i).cigar_type;
    int num = cigar_list.cigars.at(i).num;
    BamTools::CigarOp cigar_op(cigar_type, num);
    cigar_data.push_back(cigar_op);
  }
  aligned_read->cigar_ops = cigar_data;
  // Check alignment quality
  if (sw_score < min_sw_score) {
    return false;
  }
  return true;
}

ReadContainer::~ReadContainer() {
  if (output_bams) {
    delete writer_reads;
    delete writer_filtered;
  }
}
