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
  if (ID.empty())
    errx(1, "Error: found empty ID in FASTA file '%s' line %zu",
         filename.c_str(), current_line);

  if (valid_nucleotides_string(ID))
    errx(1, "Error: found Multi-lined FASTA sequence in file " \
         "'%s' line %zu (this program requires single-lined FASTA files)",
         filename.c_str(), current_line);

  if (ID.at(0) != '>')
    errx(1, "Error: found Invalid FASTA ID in file '%s' line " \
         "%zu (expected '>' character)",
         filename.c_str(), current_line);

  // If we can read the ID, but not the nucleotides,
  // This is a problematic FASTA file
  current_line++;

  if (!getline(input_stream, nuc)) {
    errx(1, "Error reading nucleotide line from FASTA file '%s' line %zu",
         filename.c_str(), current_line);
  }

  if (nuc.empty())
    errx(1, "Error: found empty nucleotide line in FASTA file '%s'" \
         "line %zu",
         filename.c_str(), current_line);
  if (!valid_nucleotides_string(nuc))
    errx(1, "Error: found invalid nucleotide line in FASTA file '%s' line %zu",
         filename.c_str(), current_line);

  read->ID = ID.substr(1);
  read->nucleotides = nuc;
  read->quality_scores = string(nuc.length(), 'N');
  read->orig_nucleotides = nuc;
  read->orig_qual = read->quality_scores;
  read->paired = false;
  return true;
}