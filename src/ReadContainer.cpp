/*
 * Author: Melissa Gymrek 2012
 */

#include <err.h>
#include <error.h>
#include <iostream>

#include "ReadContainer.h"
#include "runtime_parameters.h"

using namespace std;

ReadContainer::ReadContainer(){};

void ReadContainer::AddReadsFromFile(string bamfile) {
  if (!reader.Open(bamfile)){
    errx(1, "Could not open bam file");
  }
  std::string header_text = reader.GetHeaderText();
  user_defined_arguments = header_text;
  //  const BamTools::SamHeader header = reader.GetHeader();
  /*  user_defined_arguments = header.HasComments()?
    header.Comments.at(0):
    "# lobSTR alignment parameters not available"; // for allelotyper store header here
  */
  const BamTools::RefVector references = reader.GetReferenceData();

  BamTools::BamAlignment aln;
  while (reader.GetNextAlignment(aln)) {
    AlignedRead aligned_read;
    // get nucleotides
    aligned_read.nucleotides = aln.QueryBases;
    // get qualities
    aligned_read.qualities = aln.Qualities;
    // get strand
    aligned_read.strand = aln.IsReverseStrand();
    // get chrom
    aligned_read.chrom = references.at(aln.RefID).RefName;
    // get msStart
    if (!aln.GetTag("XS", aligned_read.msStart)) {
      aligned_read.msStart = 0;
    }
    // get msEnd
    if (!aln.GetTag("XE", aligned_read.msEnd)) {
      aligned_read.msEnd = 0;
    }
    // get read start
    aligned_read.read_start = aln.Position;
    // get cigar
    aligned_read.cigar_ops = aln.CigarData;
    // get STR seq 
    if (!aln.GetTag("XR", aligned_read.repseq)) {
      aligned_read.repseq="";
    }
    // get partial
    if (!aln.GetTag("XP",aligned_read.partial)) {
      aligned_read.partial = 0;
    }
    // get period
    aligned_read.period = aligned_read.repseq.length();
    // get diff
    if (!aln.GetTag("XD",aligned_read.diffFromRef)) {
      aligned_read.diffFromRef = 0;
    }
    if (!include_flank) { // diff is just sum of differences in cigar
	CIGAR_LIST cigar_list;
	for (vector<BamTools::CigarOp>::const_iterator it = aligned_read.cigar_ops.begin();
	       it != aligned_read.cigar_ops.end(); it++) {
	 CIGAR cig;
	 cig.num = (*it).Length;
	 cig.cigar_type = (*it).Type;
	 cigar_list.cigars.push_back(cig);
	}
	aligned_read.diffFromRef = GetSTRAllele(aligned_read, cigar_list);
    }
    if (profile) {
      cout << aligned_read.nucleotides << " " << aligned_read.diffFromRef << " "<<  aligned_read.period << endl;
    }

    if (unit) {
      if (aligned_read.diffFromRef % aligned_read.period  != 0) continue;
    }
    if (abs(aligned_read.diffFromRef) >= max_diff_ref) continue;
    // get ref copy num
    if (!aln.GetTag("XC", aligned_read.refCopyNum)) {
      aligned_read.refCopyNum = 0;
    }

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

void ReadContainer::RemovePCRDuplicates() {
  for (map<pair<string, int>, list<AlignedRead> >::iterator
	 it = aligned_str_map_.begin();
       it != aligned_str_map_.end(); it++) {
    // map of <start pos, length> -> aligned reads list
    map<pair<int, int>, list<AlignedRead> > pcr_duplicates;
    // Group into duplicates
    for (list<AlignedRead>::const_iterator
	   it2 = it->second.begin(); it2 != it->second.end(); it2++) {
      pair<int, int> key(it2->read_start, it2->nucleotides.length());
      if (pcr_duplicates.find(key) != pcr_duplicates.end()) {
	pcr_duplicates.at(key).push_back(*it2);
      } else {
	list<AlignedRead> pcr_dup_reads;
	pcr_dup_reads.push_back(*it2);
	pcr_duplicates.insert(pair< pair<int, int>, list<AlignedRead> >
			       (key, pcr_dup_reads));
      }
    }
    // Choose one rep from each group
    list<AlignedRead> reads_after_rmdup;
    for (map<pair<int, int>, list<AlignedRead> >::const_iterator
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
  return total_quality/(float)aligned_read_list.size();
}

float ReadContainer::GetScore(const string& quality_string) {
  if (quality_string.empty()) return 0.0;
  float total_quality;
  for (size_t i = 0; i < quality_string.length(); i++) {
    total_quality +=  quality_string.at(i) - 33;
  }
  return total_quality/(float)quality_string.length();
}




int ReadContainer::GetSTRAllele(const AlignedRead& aligned_read, const CIGAR_LIST& cigar_list) {
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

ReadContainer::~ReadContainer(){};
