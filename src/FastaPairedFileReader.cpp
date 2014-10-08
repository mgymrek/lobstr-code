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
#include <string>

#include "src/common.h"
#include "src/FastaPairedFileReader.h"
#include "src/runtime_parameters.h"
#include "src/ZippedFastaFileReader.h"

using namespace std;

FastaPairedFileReader::FastaPairedFileReader(const string& _filename1,
                                             const string& _filename2) {
  if (gzip) {
    _reader1 = new ZippedFastaFileReader(_filename1);
    _reader2 = new ZippedFastaFileReader(_filename2);
  } else {
    _reader1 = new FastaFileReader(_filename1);
    _reader2 = new FastaFileReader(_filename2);
  }
}

bool FastaPairedFileReader::GetNextRecord(ReadPair* read_pair) {
  read_pair->reads.clear();
  MSReadRecord read1;
  MSReadRecord read2;
  if (_reader1->GetNextRead(&read1)) {
    read1.paired = true;
    read_pair->reads.push_back(read1);
  } else {
    return false;
  }
  if (_reader2->GetNextRead(&read2)) {
    read2.paired = true;
    read_pair->reads.push_back(read2);
  } else {
    return false;
  }
  return true;
}

// Does not apply to paired read readers
bool FastaPairedFileReader::GetNextRead(MSReadRecord* read) {
  stringstream msg;
  msg << "[FastqPairedFileReader]: processing " << read->ID
      << ". This should not happen";
  PrintMessageDieOnError(msg.str(), ERROR);
  return false;
}
