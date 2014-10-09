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

#ifndef SRC_COMMON_H__
#define SRC_COMMON_H__

#include <stdio.h>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <sstream>
#include <vector>

#include "src/bwtaln.h"
#include "src/bwase.h"
#include "src/IFileReader.h"
#include "src/ReferenceSTR.h"

struct  BWT {
  bwt_t *bwt[2];
};

struct BNT {
  bntseq_t *bns;
};

struct REFSEQ {
  std::string sequence;
  std::string chrom;
  int start;
  vector<std::string> motifs;
  std::map<std::string, vector<ReferenceSTR> > ref_strs;
};

// super cool lobSTR ascii art
void PrintLobSTR();

// add option to params string
void AddOption(const std::string& optname,
               const std::string& optval,
               bool hasvalue, std::string* paramstring);

// print message to command line
// Exit if type is ERROR
enum MSGTYPE {
  ERROR = 0,
  WARNING,
  DEBUG,
  PROGRESS
};

// Output run statistics to file and to S3
void OutputRunStatistics();

// Print msg, exit if error
void PrintMessageDieOnError(const std::string& msg,
                            MSGTYPE msgtype);

// Check index version
void CheckIndexVersion();

// Debug statements
std::string GetReadDebug(const ReadPair& read_pair,
                         const std::string& detector_err,
                         const std::string& detector_msg,
                         const std::string& aln_err,
                         const std::string& aln_msg);
                         
// Get read group ID string
std::string GetReadGroup();

// trim read based on quality scores
void TrimRead(const std::string& input_nucs,
              const std::string& input_quals,
              std::string* trimmed_nucs,
              std::string* trimmed_quals,
              int cutoff);

// check if a file exists
bool fexists(const char *filename);

// check if the string contains only valid nucleotides
bool valid_nucleotides_string(const std::string &str);

// get the percentage of N's in the read
double calculate_N_percentage(const std::string& nucleotides);

// convert nucleotide to number
int nucToNumber(const char& nuc);

// check for number of repeats of most common kmer
bool CheckRepeatCount(const std::string& nucs, const size_t& k, const size_t& minlen, std::string* bestkmer);

// get the reverse complement of a nucleotide string
std::string reverseComplement(const std::string& nucs);

// get the reverse of a string
std::string reverse(const std::string& s);

// get the complement of a nucleotide
char complement(const char nucleotide);

// get the appropriate fiile reader
IFileReader* create_file_reader(const std::string& filename1,
                                const std::string& filename2);

// make sure cigar string is valid
void GenerateCorrectCigar(CIGAR_LIST* cigar_list,
                          const std::string& nucs,
                          bool* added_s,
                          bool* cigar_had_s);

// Generate the date in YYMMDD
std::string currentDateTime();

// Given number of seconds, returns a string
// describing the duration, in format "(DD days and) HH:MM:SS"
// (e.g. "3 days and 13:56:33" or just "23:00:33")
// The 'days' portion will appear only if the duration is longer than a day.
std::string GetDurationString(const size_t duration);

// Prints Running time information
void OutputRunningTimeInformation(const size_t start_time,
                                  const size_t processing_start_time,
                                  const size_t end_time,
                                  const size_t num_threads,
                                  const size_t units_processed);

// Replace string method
std::string string_replace(std::string src,
                           const std::string& target,
                           const std::string& replace);

std::vector<std::string> &split(const std::string &s, char delim,
                                std::vector<std::string> &elems);
std::string GetTime();

#endif  // SRC_COMMON_H__
