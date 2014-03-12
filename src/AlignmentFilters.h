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


#ifndef SRC_ALIGNMENTFILTERS_H
#define SRC_ALIGNMENTFILTERS_H

#include <string>

#include "src/AlignedRead.h"

using namespace std;

/* Length of perfect base matches at 5' and 3' end of read. */
pair<int,int> GetNumEndMatches(AlignedRead* aln, const string& ref_seq, int ref_seq_start);

/* Minimum distance from 5' and 3' end of reads to first indel or other end of read. */
pair<int,int> GetEndDistToIndel(AlignedRead& aln);

#endif
