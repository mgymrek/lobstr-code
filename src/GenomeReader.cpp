/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#include <err.h>
#include <error.h>
#include <iostream>
#include <stdlib.h>

#include "GenomeReader.h"
#include "runtime_parameters.h"

GenomeReader::GenomeReader(const std::string& _filename) :
  TextFileReader(_filename) {
  string line;
  string chrom;
  string sequence;
  getline(input_stream,line);
  chrom = line.substr(1);
  sequence = "";
  while (getline(input_stream,line)) {
    if (line.at(0) == '>') {
      genome.insert(pair<string,string>(chrom,sequence));
      lengths.insert(pair<string,int>(chrom,sequence.length()));
      chrom = line.substr(1);
      sequence = "";
    } else{
      sequence.append(line);
    }
  }
  if (chrom.size() <= 7) {
    genome.insert(pair<string,string>(chrom,sequence));
    lengths.insert(pair<string,int>(chrom,sequence.length()));
  }
}

bool GenomeReader::GetGenomeCoords(const string& chrom, int start,
				   int end, string* nucs) {
  if (end <=start) {cerr << "end coord must be greater than start" << endl;}
  if (genome.find(chrom) == genome.end()) {
    cerr << "chrom " << chrom << " not found in genome" << endl;
    return false;
  }
  if (end >= lengths.at(chrom)) {
    return false;
  }
  if (start < 0) {
    return false;
  }
  *nucs = genome.at(chrom).substr(start,end-start);
  return true;
}

GenomeReader::~GenomeReader() {}

bool GenomeReader::GetNextRecord(MSReadRecord* read) {
  return false;
}
