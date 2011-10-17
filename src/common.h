/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef __COMMON_FUNCIONS_H__
#define __COMMON_FUNCIONS_H__

#include <cmath>
#include <fftw3.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <istream>
#include <string>
#include <sstream>
#include <vector>

#include "bwtaln.h"
#include "bwase.h"
#include "IFileReader.h"

struct  BWT{
  bwt_t *bwt[2];
};

struct BNT {
  bntseq_t *bns;
};


// get the canonicalized form of the repeat
void getMSSeq(const std::string& nucs, int k, std::string* repeat);

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

// get the complement of a nucleotide
static char complement(const char nucleotide);

// get the canonical MS sequence
void getCanonicalMS(const std::string& msnucs, std::string* canonical);

// get the appropriate fiile reader
IFileReader* create_file_reader(const std::string& filename);

// debugging functions
std::string fftw_complex_to_string(fftw_complex v);

// print a matlab vector for debugging
template<typename TYPE>
void debug_print_matlab_vector(const std::vector<TYPE> &v, const std::string& varname) {
  std::cerr << varname << " = [ " ;
  
  if (v.size()>8)
    std::cerr << std::endl ;
  
  int col=0;
  std::streamsize old_precision =std::cerr.precision();
  std::cerr.precision(4);
  std::ios_base::fmtflags old_flags = std::cerr.flags();
  std::cerr.setf(std::ios::fixed,std::ios::floatfield);
  std::streamsize old_width = std::cerr.width();
  for ( size_t i =0 ; i<v.size(); ++i ) {
    std::cerr << std::setw(8) << v[i];
    col++;
    if (col==8) {
      col=0;
      std::cerr << " ..." << std::endl;
    } else
      std::cerr << " " ;
  }
  std::cerr << "]" << std::endl;
  
  std::cerr.precision(old_precision);
  std::cerr.flags(old_flags);
  std::cerr.width(old_width);
}

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);

#endif /* __COMMON_FUNCIONS_H__ */
