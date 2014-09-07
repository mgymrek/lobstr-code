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

#include <fftw3.h>
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

// count number of occurrences of char in string
size_t count(const std::string& s, const char& c);

// compare seqs lexicographically
std::string getFirstString(const std::string& seq1, const std::string& seq2);

// convert chromosome string to a number (e.g. "chr2" to 2)
int GetChromNumber(std::string chromosome);

// get the average quality score from a group of quality scores
float GetAverageQualityScore(const std::vector<std::string>& qualities);

// get the average quality score from single read
float GetQualityScore(const std::string& quality_score);

// check if a file exists
bool fexists(const char *filename);

// check if the string contains only valid nucleotides
bool valid_nucleotides_string(const std::string &str);

// determine which nucleotide is overrepresented
std::string OneAbundantNucleotide(const std::string& nuc, float perc_threshold);

// determine max length run of single nucleotide
int CountAbundantNucRuns(const std::string& nuc, char abundant_nuc);

// determine the magnitude of a complex number
inline double magnitude(const fftw_complex n) {
  return std::sqrt( n[0]*n[0] + n[1]*n[1] );
}

// get the percentage of N's in the read
double calculate_N_percentage(const std::string& nucleotides);

// convert nucleotide to number
int nucToNumber(const char& nuc);

// get the reverse complement of a nucleotide string
std::string reverseComplement(const std::string& nucs);

// get the reverse fo a string
std::string reverse(const std::string& s);

// get the complement of a nucleotide
char complement(const char nucleotide);

// Generate all nucleotide kmers of a certain size
void GenerateAllKmers(int size, std::vector<std::string>* kmers);

// get the minimum cyclic permutation of the provided sequence 
std::string getMinPermutation(const std::string& msnucs);

// get the canonical version of the smallest repeating subunit
std::string getCanonicalRepeat(const std::string& msnucs);

// get the appropriate fiile reader
IFileReader* create_file_reader(const std::string& filename1,
                                const std::string& filename2);

// Generate s3cmd command to get file to process
// TODO(mgymrek) eventually replace this with an
// iostream type object that can read s3 files in chunks
std::string GenerateS3Command(const std::string& bucket,
                              const std::string& filename,
                              const std::string& configfile);

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

// debugging functions
std::string fftw_complex_to_string(fftw_complex v);

// print a matlab vector for debugging
template<typename TYPE>
void debug_print_matlab_vector(const std::vector<TYPE> &v,
                               const std::string& varname) {
  std::cerr << varname << " = [ ";

  if (v.size()>8) {
    std::cerr << std::endl;
  }

  int col = 0;
  std::streamsize old_precision =std::cerr.precision();
  std::cerr.precision(4);
  std::ios_base::fmtflags old_flags = std::cerr.flags();
  std::cerr.setf(std::ios::fixed, std::ios::floatfield);
  std::streamsize old_width = std::cerr.width();
  for ( size_t i = 0 ; i < v.size(); ++i ) {
    std::cerr << std::setw(8) << v[i];
    col++;
    if (col == 8) {
      col = 0;
      std::cerr << " ..." << std::endl;
    } else {
      std::cerr << " ";
    }
  }
  std::cerr << "]" << std::endl;

  std::cerr.precision(old_precision);
  std::cerr.flags(old_flags);
  std::cerr.width(old_width);
}

std::vector<std::string> &split(const std::string &s, char delim,
                                std::vector<std::string> &elems);
std::string GetTime();

#endif  // SRC_COMMON_H__
