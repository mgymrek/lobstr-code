/*
Copyright (C) 2014 Thomas Willems <twillems@mit.edu>

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


#ifndef SRC_FILTER_COUNTER_H_
#define SRC_FILTER_COUNTER_H_

#include <stdint.h>

#include <string>

class FilterCounter {
 private:
  uint64_t* counts;
  
 public:
  const static int NUM_FILTERS     = 13;

  // Various filter types
  const static int NOT_UNIT        = 0;
  const static int DIFF_FROM_REF   = 1;
  const static int MAPPING_QUALITY = 2;
  const static int MATE_DIST       = 3;
  const static int ALLELE_SIZE     = 4;
  const static int CONTAINS_N_BASE = 5;
  const static int SPANNING_AMOUNT = 6;
  const static int NUM_END_MATCHES = 7;
  const static int NOT_MAXIMAL_END = 8;
  const static int BP_BEFORE_INDEL = 9;
  const static int CONTAINS_CLIP   = 10;
  const static int MAPQ0           = 11;
  const static int UNFILTERED      = 12;

  FilterCounter();
  
  /* Increment the count for the provided filter type */
  void increment(const int type);
  
  /* Returns the name of the filter associated with the provided filter type */
  std::string GetFilterType(const int type);

  /* Returns the count associated with the provided filter type */
  uint64_t GetFilterCount(const int type);

  ~FilterCounter();
};

#endif
