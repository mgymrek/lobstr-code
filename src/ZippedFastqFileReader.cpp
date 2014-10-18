/*
Copyright (C) 2011-2014 Melissa Gymrek <mgymrek@mit.edu>

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
#include "src/runtime_parameters.h"
#include "src/ZippedFastqFileReader.h"

using namespace std;

ZippedFastqFileReader::ZippedFastqFileReader(const std::string& _filename)
  : ZippedTextFileReader(_filename) {}

bool ZippedFastqFileReader::GetNextRecord(ReadPair* read_pair) {
  read_pair->reads.clear();
  MSReadRecord single_read;
  if (GetNextRead(&single_read)) {
    read_pair->reads.push_back(single_read);
    return true;
  } else {
    return false;
  }
}

bool ZippedFastqFileReader::GetNextRead(MSReadRecord* read) {
  string ID;
  string nuc;
  string ID2;
  string quality;

  // First line = ID
  // If no more lines, this is EOF
  current_line++;
  if (!getline(input_stream, ID))
    return false;

  // Minimal input validation
  if (ID.empty()) {
    stringstream msg;
    msg << "Found empty ID in FASTQ file " << filename << " line " << current_line;
    PrintMessageDieOnError(msg.str(), ERROR);
  }
  if (ID.at(0) != '@') {
    stringstream msg;
    msg << "Found Invalid FASTQ ID in file "
        << filename << " line " << current_line
        << " (expected '@' character)";
    PrintMessageDieOnError(msg.str(), ERROR);
  }

  // Second line = Nucleotides
  // If we can read the ID, but not the nucleotides,
  // This is a problematic FASTQ file
  current_line++;
  if (!getline(input_stream, nuc)) {
    stringstream msg;
    msg << "Problem reading nucleotide line from FASTQ file "
        << filename << " line " << current_line;
    PrintMessageDieOnError(msg.str(), ERROR);
  }

  if (nuc.empty()) {
    stringstream msg;
    msg << "Found empty nucleotide line from FASTQ file "
        << filename << " line " << current_line;
    PrintMessageDieOnError(msg.str(), ERROR);
  }
  if (!valid_nucleotides_string(nuc)) {
    stringstream msg;
    msg << "Found invalid nucleotide line from FASTQ file "
        << filename << " line " << current_line;
    PrintMessageDieOnError(msg.str(), ERROR);
  }
  // Third line = Second ID (ignored)
  current_line++;
  if (!getline(input_stream, ID2)) {
    stringstream msg;
    msg << "Problem reading second ID line from FASTQ file "
        << filename << " line " << current_line;
    PrintMessageDieOnError(msg.str(), ERROR);
  }

  // Fourth line = Quality scores (must be ASCII quality scores, not numeric)
  current_line++;
  if (!getline(input_stream, quality)) {
    stringstream msg;
    msg << "Problem reading quality scores line from FASTQ file "
        << filename << " line " << current_line;
    PrintMessageDieOnError(msg.str(), ERROR);
  }

  if (quality.length() != nuc.length()) {
    stringstream msg;
    msg << "Mismatching number of nucleotides and quality scores in FASTQ file "
        << filename << " line " << current_line;
    PrintMessageDieOnError(msg.str(), ERROR);
  }
  read->ID = ID.substr(1);
  string trim_nucs;
  string trim_qual;
  TrimRead(nuc, quality, &trim_nucs, &trim_qual, QUAL_CUTOFF);
  read->nucleotides = trim_nucs;
  read->quality_scores = trim_qual;
  read->orig_nucleotides = trim_nucs;
  read->orig_qual = trim_qual;
  read->paired = false;
  return true;
}
