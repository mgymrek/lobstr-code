/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
 */

#ifndef __CIGAR_H__
#define __CIGAR_H__

#include <string>
#include <vector>

struct CIGAR {
  int num;
  char cigar_type;
};

struct CIGAR_LIST {
  std::vector<CIGAR> cigars;
  std::string cigar_string;
};

#endif /* __CIGAR_H__ */
