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


#ifndef SRC_NWNOREFENDPENALTY_H__
#define SRC_NWNOREFENDPENALTY_H__

#include <stdlib.h>
#include <stdio.h>

#include <algorithm>
#include <string>
#include <vector>

#include "src/api/BamAux.h"

namespace NWNoRefEndPenalty { 
  void Align(const std::string& ref_seq, 
	     const std::string& read_seq,
	     std::string& ref_seq_al, 
	     std::string& read_seq_al,
	     float* score, 
	     std::vector<BamTools::CigarOp>& cigar_list);

}
#endif

