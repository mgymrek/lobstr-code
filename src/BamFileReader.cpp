/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#include <err.h>
#include <iostream>

#include "api/BamReader.h"
#include "BamFileReader.h"
#include "common.h"


using namespace std;
using BamTools::BamReader;
using BamTools::BamAlignment;

BamFileReader::BamFileReader(const std::string& _filename) {
  if (!reader.Open(_filename)) {
    errx(1, "Could not open bam file");
  }
}

bool BamFileReader::GetNextRecord(MSReadRecord* read) {
 BamTools:BamAlignment aln;
  // check if any lines left
  if (!reader.GetNextAlignment(aln)) {
    return false;
  }
  
  string ID;
  string nuc;
  string ID2;
  string quality;
  
  read->ID = aln.Name;
  read->nucleotides = aln.QueryBases;
  read->quality_scores = aln.Qualities;
  read->orig_nucleotides = aln.QueryBases;
  read->orig_qual = aln.Qualities;
  return true;
}
