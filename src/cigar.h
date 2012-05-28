/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
 */

#ifndef SRC_CIGAR_H__
#define SRC_CIGAR_H__

#include <string>
#include <sstream>
#include <vector>

struct CIGAR {
  int num;
  char cigar_type;
};

struct CIGAR_LIST {
  std::vector<CIGAR> cigars;
  std::string cigar_string;
  void ResetString() {
    std::stringstream new_cigar;
    for (std::vector<CIGAR>::const_iterator it = cigars.begin();
         it != cigars.end(); it++) {
      new_cigar << it->num << it->cigar_type;
    }
    cigar_string = new_cigar.str();
  }
};


#endif  // SRC_CIGAR_H__
