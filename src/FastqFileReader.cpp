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
#include "src/FastqFileReader.h"
#include "src/runtime_parameters.h"

using namespace std;

FastqFileReader::FastqFileReader(const std::string& _filename)
  : TextFileReader(_filename) {}

bool FastqFileReader::GetNextRecord(ReadPair* read_pair) {
  read_pair->reads.clear();
  MSReadRecord single_read;
  if (GetNextRead(&single_read)) {
    read_pair->reads.push_back(single_read);
    return true;
  } else {
    return false;
  }
}

bool FastqFileReader::GetNextRead(MSReadRecord* read) {
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
  if (ID.empty())
    errx(1, "Error: found empty ID in FASTQ file '%s' line %zu",
         filename.c_str(), current_line);

  if (ID.at(0) != '@')
    errx(1, "Error: found Invalid FASTQ ID in file '%s' "\
         "line %zu (expected '@' character)",
         filename.c_str(), current_line);

  // Second line = Nucleotides
  // If we can read the ID, but not the nucleotides,
  // This is a problematic FASTQ file
  current_line++;
  if (!getline(input_stream, nuc)) {
    errx(1, "Error reading nucleotide line from FASTQ " \
         "file '%s' line %zu",
         filename.c_str(), current_line);
  }

  if (nuc.empty())
    errx(1, "Error: found empty nucleotide line in FASTQ " \
         "file '%s' line %zu",
         filename.c_str(), current_line);
  if (!valid_nucleotides_string(nuc))
    errx(1, "Error: found invalid nucleotide line in FASTQ " \
         "file '%s' line %zu",
         filename.c_str(), current_line);

  // Third line = Second ID (ignored)
  current_line++;
  if (!getline(input_stream, ID2)) {
    errx(1, "Error reading second ID line from FASTQ file '%s' line %zu",
         filename.c_str(), current_line);
  }

  // Fourth line = Quality scores (must be ASCII quality scores, not numeric)

  current_line++;
  if (!getline(input_stream, quality)) {
    errx(1, "Error reading Quality scores line from FASTQ " \
         "file '%s' line %zu",
         filename.c_str(), current_line);
  }

  if ( quality.length() != nuc.length() )
    errx(1, "Error: Mismatching number of nucleotides (%zu) " \
         "and quality-scores (%zu) in "                      \
         "FASTQ file '%s' line '%zu'",
         nuc.length(), quality.length(), filename.c_str(), current_line);

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
