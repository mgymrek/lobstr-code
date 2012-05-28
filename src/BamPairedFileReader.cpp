/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#include <err.h>
#include "common.h"
#include "BamPairedFileReader.h"

using namespace std;

BamPairedFileReader::BamPairedFileReader(const string& _filename) {
  if (!reader.Open(_filename)) {
    errx(1, "Could not open bam file");
  }
}

bool BamPairedFileReader::GetNextRecord(ReadPair* read_pair) {
  // TODO 
  errx(1, "Paired end bam file input not yet implemented.");
  return false;
}

// Does not apply to paired read readers
bool BamPairedFileReader::GetNextRead(MSReadRecord* read) {
  return false;
}
