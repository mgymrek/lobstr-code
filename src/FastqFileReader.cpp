/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#include <err.h>
#include <error.h>
#include <iostream>

#include "common.h"
#include "FastqFileReader.h"

using namespace std;

FastqFileReader::FastqFileReader(const std::string& _filename) :
  TextFileReader(_filename) {}

bool FastqFileReader::GetNextRecord(MSReadRecord* read) {
  string ID;
  string nuc;
  string ID2;
  string quality;
  
  //
  // First line = ID
  //
  
  //If no more lines, this is EOF
  current_line++;
  if (!getline(input_stream, ID))
    return false;
  
  //Minimal input validation
  if (ID.empty())
    errx(1,"Error: found empty ID in FASTQ file '%s' line %zu",
	 filename.c_str(), current_line);
  
  if (ID.at(0) != '@')
    errx(1,"Error: found Invalid FASTQ ID in file '%s' line %zu (expected '@' character)",
	 filename.c_str(), current_line);
  
  //
  // Second line = Nucleotides
  //
  
  //If we can read the ID, but not the nucleotides,
  //This is a problematic FASTQ file
  current_line++;
  if (!getline(input_stream, nuc)) {
    errx(1,"Error reading nucleotide line from FASTQ file '%s' line %zu",
	 filename.c_str(), current_line);
  }
  
  if (nuc.empty())
    errx(1,"Error: found empty nucleotide line in FASTQ file '%s' line %zu",
	 filename.c_str(), current_line);
  if (!valid_nucleotides_string(nuc))
    errx(1,"Error: found invalid nucleotide line in FASTQ file '%s' line %zu",
	 filename.c_str(), current_line);
  
  //
  // Third line = Second ID (ignored)
  //
  current_line++;
  if (!getline(input_stream, ID2)) {
    errx(1,"Error reading second ID line from FASTQ file '%s' line %zu",
	 filename.c_str(), current_line);
  }
  
  //
  // Fourth line = Quality scores (must be ASCII quality scores, not numeric)
  //
  current_line++;
  if (!getline(input_stream, quality)) {
    errx(1,"Error reading Quality scores line from FASTQ file '%s' line %zu",
	 filename.c_str(), current_line);
  }
  
  if ( quality.length() != nuc.length() )
    errx(1,"Error: Mismatching number of nucleotides (%zu) and quality-scores (%zu) in " \
	 "FASTQ file '%s' line '%zu'",
	 nuc.length(), quality.length(), filename.c_str(), current_line);
  
  read->ID = ID.substr(1);
  read->nucleotides = nuc;
  read->quality_scores = quality;
  return true;
}
