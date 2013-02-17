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
#include "src/ZippedFastaFileReader.h"

using namespace std;

ZippedFastaFileReader::ZippedFastaFileReader(const string& _filename)
  : ZippedTextFileReader(_filename) {}

bool ZippedFastaFileReader::GetNextRecord(ReadPair* read_pair) {
  read_pair->reads.clear();
  MSReadRecord single_read;
  if (GetNextRead(&single_read)) {
    read_pair->reads.push_back(single_read);
    return true;
  } else {
    return false;
  }
}

bool ZippedFastaFileReader::GetNextRead(MSReadRecord* read) {
  string ID;
  string nuc;
  // If no more lines, this is EOF
  current_line++;
  if (!getline(input_stream, ID)) {
    return false;
  }
  // Minimal input validation
  if (ID.empty()) {
    stringstream msg;
    msg << "Found empty ID in FASTA file " << filename << " line " << current_line;
    PrintMessageDieOnError(msg.str(), ERROR);
  }
  if (valid_nucleotides_string(ID)) {
    stringstream msg;
    msg << "Found multi-lined fasta sequence in file "
        << filename << " line " << current_line
        << " (requires single-lined FASTA)";
    PrintMessageDieOnError(msg.str(), ERROR);
  }
  if (ID.at(0) != '>') {
    stringstream msg;
    msg << "Found Invalid FASTA ID in file "
        << filename << " line " << current_line
        << " (expected '>' character)";
    PrintMessageDieOnError(msg.str(), ERROR);
  }
  // If we can read the ID, but not the nucleotides,
  // This is a problematic FASTA file
  current_line++;

  if (!getline(input_stream, nuc)) {
    stringstream msg;
    msg << "Problem reading nucleotide line from FASTA file "
        << filename << " line " << current_line;
    PrintMessageDieOnError(msg.str(), ERROR);
  }

  if (nuc.empty()) {
    stringstream msg;
    msg << "Found empty nucleotide line from FASTA file "
        << filename << " line " << current_line;
    PrintMessageDieOnError(msg.str(), ERROR);
  }
  if (!valid_nucleotides_string(nuc)) {
    stringstream msg;
    msg << "Found invalid nucleotide line from FASTA file "
        << filename << " line " << current_line;
    PrintMessageDieOnError(msg.str(), ERROR);
  }
  read->ID = ID.substr(1);
  read->nucleotides = nuc;
  read->quality_scores = string(nuc.length(), 'N');
  read->orig_nucleotides = nuc;
  read->orig_qual = read->quality_scores;
  read->paired = false;
  return true;
}
