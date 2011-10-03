/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#include <boost/algorithm/string.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "common.h"
#include "AlignmentFileReader.h"
#include "runtime_parameters.h"

using namespace std;

AlignmentFileReader::AlignmentFileReader(const std::string& _filename) :
  TextFileReader(_filename) {
  // read first line
  string line;
  current_line++;
  getline(input_stream, line);
}

bool AlignmentFileReader::GetNextRecord(MSReadRecord* rec) {
  // declare variables to put in msrecord
  string line;
  // if no more lines, this is EOF
  current_line++;
  if(!getline(input_stream,line))
    return false;
  // populate the msrecord
  vector<string> items;
  split(line, '\t', items);
  rec->ID = items[0];
  rec->nucleotides = items[1];
  rec->quality_scores = items[2];
  rec->ms_start = atoi(items[3].c_str());
  rec->ms_end = atoi(items[4].c_str());
  rec->ms_repeat_best_period = atoi(items[5].c_str());
  rec->left_flank_nuc = items[6];
  rec->detected_ms_region_nuc = items[7];
  rec->right_flank_nuc = items[8];
  rec->chrom = items[9];
  rec->msStart = atoi(items[10].c_str());
  rec->msEnd = atoi(items[11].c_str());
  rec->msRepeat = items[12];
  rec->refCopyNum = atof(items[13].c_str());
  rec->lStart = atoi(items[14].c_str());
  rec->lEnd = atoi(items[15].c_str());
  rec->rStart = atoi(items[16].c_str());
  rec->rEnd = atoi(items[17].c_str());
  rec->diffFromRef = atof(items[18].c_str());
  rec->reverse = (items[19] == "1");
  if (items.size() > 21) rec->name = items[21];
  return true;
}


AlignmentFileReader::~AlignmentFileReader() {}
