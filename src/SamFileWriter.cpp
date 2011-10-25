/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#include "common.h"
#include "SamFileWriter.h"

using BamTools::BamWriter;
using BamTools::RefData;
using BamTools::RefVector;
using BamTools::BamAlignment;

SamFileWriter::SamFileWriter(const std::string& _filename,
			     const map<string, int>& _chrom_sizes) {//:TextFileWriter(_filename) {
  chrom_sizes = _chrom_sizes;
  // open bam writer
  const string header = "";
  RefVector ref_vector;
  for (map<string, int>::const_iterator it = chrom_sizes.begin();
       it != chrom_sizes.end(); ++it) {
    const string& name = it->first;
    int length = it->second;
    RefData ref_data(length);
    ref_data.RefName = name;
    ref_vector.push_back(ref_data);
  }
  if (!writer.Open(_filename, header, ref_vector)) {
    cerr << "Could not open bam file" << endl;
  }
}

void SamFileWriter::WriteRecord(const MSReadRecord &msread) {
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
  bam_alignment.SetIsMateUnmapped(false);
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
    BamTools::CigarOp cigar_op_l('M', msread.msStart - lpos + 1);
    BamTools::CigarOp cigar_op_m('I', msread.diffFromRef);
    BamTools::CigarOp cigar_op_r('M', (read.length() - msread.diffFromRef - (msread.msStart-lpos+1)));
    cigar_data.push_back(cigar_op_l);
    cigar_data.push_back(cigar_op_m);
    cigar_data.push_back(cigar_op_r);
  } else {
    BamTools::CigarOp cigar_op_l('M', (msread.msStart-lpos +1));
    BamTools::CigarOp cigar_op_m('D', (-1*msread.diffFromRef));
    BamTools::CigarOp cigar_op_r('M', (read.length() - (msread.msStart-lpos+1)));
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
  writer.SaveAlignment(bam_alignment);
}

SamFileWriter::~SamFileWriter() {
  writer.Close();
}


/*
  // set the params
  qual = msread.quality_scores;
  qname = msread.ID;
  flag = 16;
  rname = msread.chrom;
  lpos = msread.lStart;
  rpos = msread.rStart;
  mapq = 255;
  mrnm = "*";
  mpos = 0;
  isize = 0;
  // write output
  if (msread.diffFromRef == 0) {
    output_stream <<
      qname << "\t" <<
      flag << "\t" <<
      rname << "\t" <<
      lpos << "\t" <<
      mapq << "\t" <<
      (lseq.length() + middle.length()+rseq.length()) << "M" << "\t" <<
      mrnm << "\t"<<
      mpos << "\t" <<
      isize << "\t" <<
      read << "\t" <<
      qual;
  } else if (msread.diffFromRef > 0) {
    output_stream <<
      qname << "\t" <<
      flag << "\t" <<
      rname << "\t" <<
      lpos << "\t" <<
      mapq << "\t" <<
      (msread.msStart-lpos +1) << "M" << msread.diffFromRef << "I" << (read.length() - msread.diffFromRef - (msread.msStart-lpos+1)) << "M" <<"\t" <<
      mrnm << "\t"<<
      mpos << "\t" <<
      isize << "\t" <<
      read << "\t" <<
      qual;
  } else{
    output_stream <<
      qname << "\t" <<
      flag << "\t" <<
      rname << "\t" <<
      lpos << "\t" <<
      mapq << "\t" <<
      (msread.msStart-lpos +1) << "M" << (-1*msread.diffFromRef) << "D" << (read.length() - (msread.msStart-lpos+1)) << "M" <<"\t" <<
      mrnm << "\t"<<
      mpos << "\t" <<
      isize << "\t" <<
      read << "\t" <<
      qual;
  }  
  output_stream << "\t"
		<< "XP:Z:" << msread.chrom << ":" << msread.msStart << "-" << msread.msEnd
		<< "\t"
		<< "XL:f:" << str_length << endl;
  float str_length = msread.refCopyNum +
    float(msread.diffFromRef)/float(msread.msRepeat.length());

*/
