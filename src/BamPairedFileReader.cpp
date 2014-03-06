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
#include <map>
#include <utility>

#include "src/common.h"
#include "src/BamPairedFileReader.h"
#include "src/runtime_parameters.h"

using BamTools::BamReader;
using BamTools::BamAlignment;

using namespace std;

BamPairedFileReader::BamPairedFileReader(const std::string& _filename) {
  if (!reader.Open(_filename)) {
    PrintMessageDieOnError("Could not open bam file", ERROR);
  }
}


bool BamPairedFileReader::GetNextRecord(ReadPair* read_pair) {
  read_pair->reads.clear();
  MSReadRecord single_read;
  MSReadRecord mate;
  int64_t bam_file_position;
  if (GetNextRead(&single_read)) {
    read_pair->reads.push_back(single_read);
    if (!single_read.paired) {
      return true;
    } else {
      if (GetNextReadMate(&mate, &bam_file_position)) {
        if (mate.ID == single_read.ID) {
          read_pair->reads.push_back(mate);
          return true;
        } else {
          PrintMessageDieOnError("Could not find pair for " + single_read.ID +
                                 ". Is the bam file sorted by read name?", WARNING);
          // set single read to not paired
          read_pair->reads.at(0).paired = false;
          // back up by one read
          if (!reader.Seek(bam_file_position)) {
            PrintMessageDieOnError("Could not rewind bam file reader", ERROR);
          }
          return true;
        }
      }
    }
  } else {
    return false;
  }
  return false;
}

bool BamPairedFileReader::GetNextReadMate(MSReadRecord* read,
                                          int64_t* bam_file_position) {
  // Get position in case we need to rewind
  *bam_file_position = reader.Tell();
  BamAlignment aln;
  // check if any lines left
  if (!reader.GetNextAlignment(aln)) {
    return false;
  }
  read->ID = aln.Name;
  // strip /1 or /2 for pairs
  if (read->ID.length() > 2) {
    int pos1 = read->ID.find("/1");
    int pos2 = read->ID.find("/2");
    if (pos1 != -1) {
      read->ID = read->ID.substr(0, pos1);
    } else if (pos2 != -1) {
      read->ID = read->ID.substr(0, pos2);
    }
  }
  string nucs = aln.QueryBases;
  string qual = aln.Qualities;
  nucs = reverseComplement(nucs);
  qual = reverse(qual);
  string trim_nucs;
  string trim_qual;
  TrimRead(nucs, qual, &trim_nucs, &trim_qual, QUAL_CUTOFF);
  read->nucleotides = trim_nucs;
  read->quality_scores = trim_qual;
  read->orig_nucleotides = trim_nucs;
  read->orig_qual = trim_qual;
  read->paired = aln.IsPaired();
  return true;
}

bool BamPairedFileReader::GetNextRead(MSReadRecord* read) {
  BamAlignment aln;
  // check if any lines left
  if (!reader.GetNextAlignment(aln)) {
    return false;
  }
  read->ID = aln.Name;
  // strip /1 or /2 for pairs
  if (read->ID.length() > 2) {
    int pos1 = read->ID.find("/1");
    int pos2 = read->ID.find("/2");
    if (pos1 != -1) {
      read->ID = read->ID.substr(0, pos1);
    } else if (pos2 != -1) {
      read->ID = read->ID.substr(0, pos2);
    }
  }
  string nucs = aln.QueryBases;
  string qual = aln.Qualities;
  string trim_nucs;
  string trim_qual;
  TrimRead(nucs, qual, &trim_nucs, &trim_qual, QUAL_CUTOFF);
  read->nucleotides = trim_nucs;
  read->quality_scores = trim_qual;
  read->orig_nucleotides = trim_nucs;
  read->orig_qual = trim_qual;
  read->paired = aln.IsPaired();
  return true;
}
