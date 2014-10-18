/*
Copyright (C) 2011-2014 Melissa Gymrek <mgymrek@mit.edu>

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

Modified from http://www.rolfmuertter.com/code/nw.h.html
Rolf Muertter,  rolf@dslextreme.com
9/5/2006
*/

#ifndef SRC_NW_H__
#define SRC_NW_H__

#include <stdlib.h>
#include <stdio.h>

#include <algorithm>
#include <string>
#include <vector>

#include "src/cigar.h"

int nw(const std::string& seq_1, const std::string& seq_2,
       std::string& seq_1_al, std::string& seq_2_al,
       int* score , CIGAR_LIST* cigar_list);

int nw_align(std::vector<int>* F, std::vector<char>* traceback,
             const std::string& seq_1, const std::string& seq_2,
             std::string& seq_1_al, std::string& seq_2_al,
             int d, int* score);

/* sw alignment with affine gap penalty */
int nw_align_ag(std::vector<int>* M, std::vector<int>* I,
                std::vector<char>* tracebackM,
                std::vector<char>* tracebackI,
                const std::string& seq_1, const std::string& seq_2,
                std::string& seq_1_al, std::string& seq_2_al,
                int* score, CIGAR_LIST* cigar_list);

void  dpm_init(std::vector<int>* F, std::vector<char>* traceback,
               int L1, int L2, int d);

int max(int f1, int f2, int f3, char * ptr);
int maxM(int f1, int f2, char * ptr);
int maxI(int f1, int f2, int f3, int f4, char * ptr);

#endif  // SRC_NW_H__

