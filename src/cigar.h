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
