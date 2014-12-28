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

#include "src/common.h"
#include "src/FilterCounter.h"

FilterCounter::FilterCounter(){
  counts = new uint64_t[NUM_FILTERS];
  for (int i = 0; i < NUM_FILTERS; i++)
    counts[i] = 0;
}
  
void FilterCounter::increment(const int type){
  if (type > NUM_FILTERS || type < 0)
    PrintMessageDieOnError("Invalid filter type", ERROR);
  counts[type]++;
}
 
std::string FilterCounter::GetFilterType(const int type){
  switch(type) {
  case NOT_UNIT:
    return "NOT_UNIT";
  case DIFF_FROM_REF:
    return "DIFF_FROM_REF";
  case MAPPING_QUALITY:
    return "MAPPING_QUALITY";
  case MATE_DIST:
    return "MATE_DIST";
  case ALLELE_SIZE:
    return "ALLELE_SIZE";
  case CONTAINS_N_BASE:
    return "CONTAINS_N_BASE";
  case SPANNING_AMOUNT:
    return "SPANNING_AMOUNT";
  case NUM_END_MATCHES:
    return "NUM_END_MATCHES";
  case NOT_MAXIMAL_END:
    return "NOT_MAXIMAL_END";
  case BP_BEFORE_INDEL:
    return "BP_BEFORE_INDEL";
  case UNFILTERED:
    return "UNFILTERED";
  default:
    PrintMessageDieOnError("Invalid filter type", ERROR);
    return "ERROR"; // this should never be reached
  }
}

uint64_t FilterCounter::GetFilterCount(const int type){
  if (type > NUM_FILTERS || type < 0)
    PrintMessageDieOnError("Invalid filter type", ERROR);
  return counts[type];
}

 
FilterCounter::~FilterCounter(){
  delete [] counts;
}


  
  
