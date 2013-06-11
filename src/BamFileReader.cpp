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

#include "api/BamReader.h"
#include "src/BamFileReader.h"
#include "src/common.h"
#include "src/runtime_parameters.h"

using BamTools::BamReader;
using BamTools::BamAlignment;

using namespace std;
// Hello Melissa
BamFileReader::BamFileReader(const std::string& _filename) {
  if (!reader.Open(_filename)) {
    PrintMessageDieOnError("Could not open bam file", ERROR);
  }
}

bool BamFileReader::GetNextRecord(ReadPair* read_pair) {
  read_pair->reads.clear();
  MSReadRecord single_read;
  if (GetNextRead(&single_read)) {
    read_pair->reads.push_back(single_read);
    return true;
  } else {
    return false;
  }
}

bool BamFileReader::GetNextRead(MSReadRecord* read) {
  BamAlignment aln;
  // check if any lines left
  if (!reader.GetNextAlignment(aln)) {
    return false;
  }
  read->ID = aln.Name;
  string trim_nucs;
  string trim_qual;
  TrimRead(aln.QueryBases, aln.Qualities, &trim_nucs, &trim_qual, QUAL_CUTOFF);
  read->nucleotides = trim_nucs;
  read->quality_scores = trim_qual;
  read->orig_nucleotides = trim_nucs;
  read->orig_qual = trim_qual;
  read->paired = false;
  return true;
}
