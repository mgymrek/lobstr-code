/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#include <err.h>
#include "common.h"
#include "FastqPairedFileReader.h"
#include "runtime_parameters.h"
#include "ZippedFastqFileReader.h"

using namespace std;

FastqPairedFileReader::FastqPairedFileReader(const string& _filename1,
					     const string& _filename2) {
  if (gzip) {
    _reader1 = new ZippedFastqFileReader(_filename1);
    _reader2 = new ZippedFastqFileReader(_filename2);
  } else {
    _reader1 = new FastqFileReader(_filename1);
    _reader2 = new FastqFileReader(_filename2);
  }
}

bool FastqPairedFileReader::GetNextRecord(ReadPair* read_pair) {
  read_pair->reads.clear();
  MSReadRecord read1;
  MSReadRecord read2;
  if (_reader1->GetNextRead(&read1)) {
    read_pair->reads.push_back(read1);
  } else {
    return false;
  }
  if (_reader2->GetNextRead(&read2)) {
    read_pair->reads.push_back(read2);
  } else {
    return false;
  }
  return true;
}

// Does not apply to paired read readers
bool FastqPairedFileReader::GetNextRead(MSReadRecord* read) {
  return false;
}
