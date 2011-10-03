/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#include <boost/algorithm/string.hpp>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "common.h"
#include "MSTableFileReader.h"
#include "runtime_parameters.h"

using namespace std;

MSTableFileReader::MSTableFileReader(const std::string& _filename,
				     const std::string& _genome_file) :
  TextFileReader(_filename) {
  gReader = new GenomeReader(_genome_file);
}

void MSTableFileReader::GetChromSizes(map<string, int>* chrom_sizes) {
  *chrom_sizes = gReader->lengths;
}

bool MSTableFileReader::GetNextRecord(MSRecord* rec) {
  // declare variables to put in msrecord
  string line;
  // if no more lines, this is EOF
  current_line++;
  if(!getline(input_stream,line))
    return false;
  // populate the msrecord
  vector<string> items;
  split(line, '\t', items);
  //boost::split(items, line, boost::is_any_of("\t"));
  rec->chrom = items.at(0);
  rec->start = atoi(items.at(1).c_str());
  rec->end = atoi(items.at(2).c_str());
  rec->extend = extend;
  string repeat = items[15];
  getCanonicalMS(repeat, &rec->repeat);
  rec->copynum = atof(items.at(5).c_str());
  if (items.size() >=20 ) {
    rec->name = items.at(19);
  } else { rec->name = "";}
  int start_left = max(rec->start-extend, 0);
  int end_right = min(rec->end + extend, gReader->lengths.at(rec->chrom));
  string left_flank;
  string right_flank;

  if (!gReader->GetGenomeCoords(rec->chrom, start_left, rec->end, &left_flank)) {
    left_flank = "";
    //    return false;
  }
  boost::to_upper(left_flank);
  rec->leftFlank = left_flank;
  if (!gReader->GetGenomeCoords(rec->chrom, rec->start, end_right, &right_flank)) {
    right_flank = "";
    //return false;
    }
  /*
  if (!gReader->GetGenomeCoords(rec->chrom, start_left, rec->start + 50 , &left_flank)) {
    left_flank = "";
    //    return false;
  }
  boost::to_upper(left_flank);
  rec->leftFlank = left_flank;
  if (!gReader->GetGenomeCoords(rec->chrom, rec->end - 50, end_right, &right_flank)) {
    right_flank = "";
    //return false;
    }*/
  boost::to_upper(right_flank);
  rec->rightFlank = right_flank;
  return true;
}

bool MSTableFileReader::GetNextRecord(MSReadRecord* msrec) { return false;}

MSTableFileReader::~MSTableFileReader() {
  delete gReader;
}
