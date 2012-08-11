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

#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <sstream>
#include <vector>

#include "src/bwtaln.h"
#include "src/bwase.h"
#include "src/IFileReader.h"

struct  BWT {
  bwt_t *bwt[2];
};

struct BNT {
  bntseq_t *bns;
};

struct REFSEQ {
  std::string sequence;
  int start;
};

// add option to params string
void AddOption(const std::string& optname,
               const std::string& optval,
               bool hasvalue, std::string* paramstring);

// trim read based on quality scores
void TrimRead(const std::string& input_nucs,
              const std::string& input_quals,
              std::string* trimmed_nucs,
              std::string* trimmed_quals,
              int cutoff);

// count number of occurrences of char in string
size_t count(const std::string& s, const char& c);

// get the canonicalized form of the repeat
bool getMSSeq(const std::string& nucs, int k, std::string* repeat);

// compare seqs lexicographically
std::string getFirstString(const std::string& seq1, const std::string& seq2);

// determine if sequence is perfect repeat
bool IsPerfectRepeat(const std::string& sequence, const std::string& repeat);

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
char OneAbundantNucleotide(const std::string& nuc, float perc_threshold);

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

// get the canonical MS sequence
void getCanonicalMS(const std::string& msnucs, std::string* canonical);

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


#endif  // SRC_COMMON_H__
