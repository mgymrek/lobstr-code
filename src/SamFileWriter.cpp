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
			     const map<string, int>& _chrom_sizes) {//:TextFileWriter(_filename) {
  chrom_sizes = _chrom_sizes;
  std::string header="@CO\t";
  header += user_defined_arguments;
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

void SamFileWriter::WriteAdjustedRecord(const MSReadRecord &msread) {
  BamAlignment bam_alignment;
  lpos = msread.read_start;
  bam_alignment.Name = msread.ID;
  bam_alignment.SetIsPaired(false);
  bam_alignment.SetIsDuplicate(false);
  bam_alignment.SetIsFailedQC(false);
  bam_alignment.SetIsFirstMate(false);
  bam_alignment.SetIsMateMapped(false);
  bam_alignment.SetIsMateReverseStrand(false);
  bam_alignment.SetIsProperPair(false);
  bam_alignment.SetIsSecondMate(false);
  bam_alignment.SetIsMapped(true);
  bam_alignment.SetIsMateMapped(true);
  bam_alignment.SetIsPrimaryAlignment(true);
  if (msread.reverse) {
    bam_alignment.SetIsReverseStrand(true);
  } else {
    bam_alignment.SetIsReverseStrand(false);
  }
  bam_alignment.Position = msread.read_start;
  bam_alignment.MapQuality = 255;
  bam_alignment.Qualities = msread.quality_scores;
  if (msread.reverse) {
    bam_alignment.QueryBases = reverseComplement(msread.nucleotides);
  } else {
    bam_alignment.QueryBases = msread.nucleotides;
  }
  // rname (ref id)
  int ref_id;
  size_t i = 0;
  for (map<string, int>::const_iterator it = chrom_sizes.begin();
       it != chrom_sizes.end(); ++it) {
    if (it->first == msread.chrom) {
      ref_id = i;
    }
    ++i;
  }
  bam_alignment.RefID = ref_id;
  // cigar
  vector<BamTools::CigarOp> cigar_data;
  for (i = 0; i < msread.cigar.size(); i++) {
    char cigar_type = msread.cigar.at(i).cigar_type;
    int num = msread.cigar.at(i).num;
    BamTools::CigarOp cigar_op(cigar_type,num);
    cigar_data.push_back(cigar_op);
  }
  bam_alignment.CigarData = cigar_data;

  // write user flags giving repeat information
  // XS: start pos of matching STR
  bam_alignment.AddTag("XS", "i", msread.msStart);
  // XE: end pos of matching STR
  bam_alignment.AddTag("XE", "i", msread.msEnd);
  // XR: STR repeat
  bam_alignment.AddTag("XR", "Z", msread.msRepeat);
  // XD: nuc diff compared to ref
  bam_alignment.AddTag("XD", "i", msread.diffFromRef);
  // XC: ref copy number
  bam_alignment.AddTag("XC", "f", msread.refCopyNum);
  // XG: repeat region
  bam_alignment.AddTag("XG","Z",msread.detected_ms_nuc);
  // XW: sw score
  bam_alignment.AddTag("XW", "i", msread.sw_score);
  // XP: partial alignment
  bam_alignment.AddTag("XP", "i", (int)msread.partial);
  // XN: name of STR repeat
  if (!msread.name.empty()) {
    bam_alignment.AddTag("XN", "Z", msread.name);
  }
  if ((msread.msStart-lpos +1) <= ((int)msread.nucleotides.length())) {
    writer.SaveAlignment(bam_alignment);
  }
}


void SamFileWriter::WriteRecord(const MSReadRecord &msread) {
  if (adjust) {
    WriteAdjustedRecord(msread);
    return;
  }
  BamAlignment bam_alignment;
  if (msread.reverse) {
    lpos = msread.rStart;
  } else {
    lpos = msread.lStart;
  }
  bam_alignment.Name = msread.ID;
  bam_alignment.SetIsPaired(false);
  bam_alignment.SetIsDuplicate(false);
  bam_alignment.SetIsFailedQC(false);
  bam_alignment.SetIsFirstMate(false);
  bam_alignment.SetIsMateMapped(false);
  bam_alignment.SetIsMateReverseStrand(false);
  bam_alignment.SetIsProperPair(false);
  bam_alignment.SetIsSecondMate(false);
  bam_alignment.SetIsMapped(true);
  bam_alignment.SetIsMateMapped(true);
  bam_alignment.SetIsPrimaryAlignment(true);
  if (msread.reverse) {
    bam_alignment.SetIsReverseStrand(true);
  }
  bam_alignment.Position = msread.lStart;
  bam_alignment.MapQuality = 255;

  if (msread.reverse) {
    lseq = reverseComplement(msread.right_flank_nuc);
    middle = reverseComplement(msread.detected_ms_region_nuc);
    rseq = reverseComplement(msread.left_flank_nuc);
  } else {
    lseq = msread.left_flank_nuc;
    middle = msread.detected_ms_region_nuc;
    rseq = msread.right_flank_nuc;
  }
  read = lseq + middle + rseq;
  bam_alignment.Qualities = msread.quality_scores;
  bam_alignment.QueryBases = read;
  // rname (ref id)
  int ref_id;
  int i = 0;
  for (map<string, int>::const_iterator it = chrom_sizes.begin();
       it != chrom_sizes.end(); ++it) {
    if (it->first == msread.chrom) {
      ref_id = i;
    }
    ++i;
  }
  bam_alignment.RefID = ref_id;
  // cigar
  vector<BamTools::CigarOp> cigar_data;
  if (msread.diffFromRef == 0) {
    BamTools::CigarOp cigar_op('M', lseq.length() + middle.length()
			       + rseq.length());
    cigar_data.push_back(cigar_op);
  } else if (msread.diffFromRef > 0) {
    size_t m1 = msread.msStart - lpos + 1;
    int ins = msread.diffFromRef;
    size_t m2 = (read.length() - msread.diffFromRef - (msread.msStart-lpos+1));
    if ((m1 + ins + m2 != read.length()) || (m1 < 0) || (m2 < 0))
      return;
    BamTools::CigarOp cigar_op_l('M', m1);
    BamTools::CigarOp cigar_op_m('I', ins);
    BamTools::CigarOp cigar_op_r('M', m2);
    cigar_data.push_back(cigar_op_l);
    cigar_data.push_back(cigar_op_m);
    cigar_data.push_back(cigar_op_r);
  } else {
    size_t m1 = (msread.msStart-lpos +1);
    int del = (-1*msread.diffFromRef);
    size_t m2 = (read.length() - (msread.msStart-lpos+1));
    if ((m1+m2 != read.length()) ||( m1 < 0) || (m2 < 0))
      return;
    BamTools::CigarOp cigar_op_l('M', m1);
    BamTools::CigarOp cigar_op_m('D', del);
    BamTools::CigarOp cigar_op_r('M', m2);
    cigar_data.push_back(cigar_op_l);
    cigar_data.push_back(cigar_op_m);
    cigar_data.push_back(cigar_op_r);
  }
  bam_alignment.CigarData = cigar_data;

  // write user flags giving repeat information
  // XS: start pos of matching STR
  bam_alignment.AddTag("XS", "i", msread.msStart);
  // XE: end pos of matching STR
  bam_alignment.AddTag("XE", "i", msread.msEnd);
  // XR: STR repeat
  bam_alignment.AddTag("XR", "Z", msread.msRepeat);
  // XD: nuc diff compared to ref
  bam_alignment.AddTag("XD", "i", msread.diffFromRef);
  // XC: ref copy number
  bam_alignment.AddTag("XC", "f", msread.refCopyNum);
  // XA: mismatch in left flank
  //bam_alignment.AddTag("XA", "i", msread.lDist);
  // XB: mismatch in right flank
  //bam_alignment.AddTag("XB", "i", msread.rDist);
  // XN: name of STR repeat
  if (!msread.name.empty()) {
    bam_alignment.AddTag("XN", "Z", msread.name);
  }
  if ((msread.msStart-lpos +1) <= ((int)msread.nucleotides.length())) {
    writer.SaveAlignment(bam_alignment);
  }
}

SamFileWriter::~SamFileWriter() {
  writer.Close();
}

