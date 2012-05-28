/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#include "common.h"
#include "runtime_parameters.h"
#include "SamFileWriter.h"

using namespace std;
using BamTools::BamWriter;
using BamTools::RefData;
using BamTools::RefVector;
using BamTools::BamAlignment;

SamFileWriter::SamFileWriter(const std::string& _filename,
			     const map<string, int>& _chrom_sizes) {
  chrom_sizes = _chrom_sizes;
  std::string header="@CO\t";
  header += user_defined_arguments;
  header +="\n";
  RefVector ref_vector;
  for (map<string, int>::const_iterator it = chrom_sizes.begin();
       it != chrom_sizes.end(); ++it) {
    const string& name = it->first;
    int length = it->second;
    RefData ref_data(name, length);
    ref_data.RefName = name;
    ref_vector.push_back(ref_data);
  }
  if (!writer.Open(_filename, header, ref_vector)) {
    cerr << "Could not open bam file" << endl;
  }
}

void SamFileWriter::WriteRecord(const ReadPair& read_pair) {
  const int& aligned_read_num = read_pair.aligned_read_num;
  const int& paired_dist = read_pair.treat_as_paired ?
    abs(read_pair.reads.at(aligned_read_num).read_start-
	read_pair.reads.at(1-aligned_read_num).read_start) : -1;
  // Write aligned read
  BamAlignment bam_alignment;
  bam_alignment.Name = read_pair.reads.at(aligned_read_num).ID;
  bam_alignment.SetIsPaired(read_pair.treat_as_paired);
  bam_alignment.SetIsDuplicate(false);
  bam_alignment.SetIsFailedQC(false);
  bam_alignment.SetIsFirstMate(read_pair.treat_as_paired);
  bam_alignment.SetIsMateMapped(read_pair.treat_as_paired);
  bam_alignment.SetIsMateReverseStrand(!read_pair.reads.at(aligned_read_num).reverse);
  bam_alignment.SetIsProperPair(read_pair.treat_as_paired);
  bam_alignment.SetIsSecondMate(false);
  bam_alignment.SetIsMapped(true);
  bam_alignment.SetIsMateMapped(read_pair.treat_as_paired);
  bam_alignment.SetIsPrimaryAlignment(true);
  bam_alignment.SetIsReverseStrand(read_pair.reads.at(aligned_read_num).reverse);
  bam_alignment.Position = read_pair.reads.at(aligned_read_num).read_start;
  bam_alignment.MapQuality = 255; // TODO change
  bam_alignment.Qualities = read_pair.reads.at(aligned_read_num).quality_scores;
  if (read_pair.reads.at(aligned_read_num).reverse) {
    bam_alignment.QueryBases = reverseComplement(read_pair.reads.at(aligned_read_num).nucleotides);
  } else {
    bam_alignment.QueryBases = read_pair.reads.at(aligned_read_num).nucleotides;
  }
  // rname (ref id)
  int ref_id;
  size_t i = 0;
  for (map<string, int>::const_iterator it = chrom_sizes.begin();
       it != chrom_sizes.end(); ++it) {
    if (it->first == read_pair.reads.at(aligned_read_num).chrom) {
      ref_id = i;
    }
    ++i;
  }
  bam_alignment.RefID = ref_id;
  // cigar
  vector<BamTools::CigarOp> cigar_data;
  for (i = 0; i < read_pair.reads.at(aligned_read_num).cigar.size(); i++) {
    char cigar_type = read_pair.reads.at(aligned_read_num).cigar.at(i).cigar_type;
    int num = read_pair.reads.at(aligned_read_num).cigar.at(i).num;
    BamTools::CigarOp cigar_op(cigar_type,num);
    cigar_data.push_back(cigar_op);
  }
  bam_alignment.CigarData = cigar_data;

  // write user flags giving repeat information
  // XS: start pos of matching STR
  bam_alignment.AddTag("XS", "i", read_pair.reads.at(aligned_read_num).msStart);
  // XE: end pos of matching STR
  bam_alignment.AddTag("XE", "i", read_pair.reads.at(aligned_read_num).msEnd);
  // XR: STR repeat
  bam_alignment.AddTag("XR", "Z", read_pair.reads.at(aligned_read_num).msRepeat);
  // XD: nuc diff compared to ref
  bam_alignment.AddTag("XD", "i", read_pair.reads.at(aligned_read_num).diffFromRef);
  // XC: ref copy number
  bam_alignment.AddTag("XC", "f", read_pair.reads.at(aligned_read_num).refCopyNum);
  // XG: repeat region
  bam_alignment.AddTag("XG","Z",read_pair.reads.at(aligned_read_num).detected_ms_nuc);
  // XW: sw score
  bam_alignment.AddTag("XW", "i", read_pair.reads.at(aligned_read_num).sw_score);
  // XP: partial alignment
  bam_alignment.AddTag("XP", "i", (int)read_pair.reads.at(aligned_read_num).partial);
  // XX: stitched
  bam_alignment.AddTag("XX", "i", (int)(!read_pair.treat_as_paired && paired));
  // XM: distance between read and mate start pos
  bam_alignment.AddTag("XM", "i", paired_dist);
  // XN: name of STR repeat
  if (!read_pair.reads.at(aligned_read_num).name.empty()) {
    bam_alignment.AddTag("XN", "Z", read_pair.reads.at(aligned_read_num).name);
  }

  writer.SaveAlignment(bam_alignment);

  // Write mate pair
  if (read_pair.treat_as_paired) {
      BamAlignment mate_alignment;
      mate_alignment.Name = read_pair.reads.at(1-aligned_read_num).ID;
      mate_alignment.SetIsPaired(true);
      mate_alignment.SetIsDuplicate(false);
      mate_alignment.SetIsFailedQC(false);
      mate_alignment.SetIsFirstMate(false);
      mate_alignment.SetIsMateMapped(true);
      mate_alignment.SetIsMateReverseStrand(read_pair.reads.at(aligned_read_num).reverse);
      mate_alignment.SetIsProperPair(true);
      mate_alignment.SetIsSecondMate(true);
      mate_alignment.SetIsMapped(true);
      mate_alignment.SetIsMateMapped(true);
      mate_alignment.SetIsPrimaryAlignment(true);
      mate_alignment.SetIsReverseStrand(read_pair.reads.at(1-aligned_read_num).reverse);
      mate_alignment.Position = read_pair.reads.at(1-aligned_read_num).read_start;
      mate_alignment.MapQuality = 255; // TODO change
      mate_alignment.Qualities = read_pair.reads.at(1-aligned_read_num).quality_scores;
      if (read_pair.reads.at(1-aligned_read_num).reverse) {
	mate_alignment.QueryBases = reverseComplement(read_pair.reads.at(1-aligned_read_num).nucleotides);
      } else {
	mate_alignment.QueryBases = read_pair.reads.at(1-aligned_read_num).nucleotides;
      }
      mate_alignment.RefID = ref_id;
      // cigar
      vector<BamTools::CigarOp> cigar_data;
      for (i = 0; i < read_pair.reads.at(1-aligned_read_num).cigar.size(); i++) {
	char cigar_type = read_pair.reads.at(1-aligned_read_num).cigar.at(i).cigar_type;
	int num = read_pair.reads.at(1-aligned_read_num).cigar.at(i).num;
	BamTools::CigarOp cigar_op(cigar_type,num);
	cigar_data.push_back(cigar_op);
      }
      mate_alignment.CigarData = cigar_data;

      // write user flags giving repeat information
      // XS: start pos of matching STR
      bam_alignment.AddTag("XS", "i", read_pair.reads.at(aligned_read_num).msStart);
      // XE: end pos of matching STR
      bam_alignment.AddTag("XE", "i", read_pair.reads.at(aligned_read_num).msEnd);
      // XR: STR repeat
      bam_alignment.AddTag("XR", "Z", read_pair.reads.at(aligned_read_num).msRepeat);
      // XC: ref copy number
      bam_alignment.AddTag("XC", "f", read_pair.reads.at(aligned_read_num).refCopyNum);
      // XN: name of STR repeat
      if (!read_pair.reads.at(aligned_read_num).name.empty()) {
	mate_alignment.AddTag("XN", "Z", read_pair.reads.at(aligned_read_num).name);
      }
      writer.SaveAlignment(mate_alignment);
  }
}

SamFileWriter::~SamFileWriter() {
  writer.Close();
}

