/*
 * nw.h for program nw.
 * taken from http://www.rolfmuertter.com/code/nw.h.html
 *   Rolf Muertter,  rolf@dslextreme.com
 *   9/5/2006
 *
 */

#include <iostream>
#include <string>
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <vector>
#include "cigar.h"


int nw(const std::string& seq_1, const std::string& seq_2, 
       std::string& seq_1_al, std::string& seq_2_al,
       bool prm, int* score , CIGAR_LIST* cigar_list);

int nw_align(std::vector<int>* F, std::vector<char>* traceback,
	     const std::string& seq_1, const std::string& seq_2, 
	     std::string& seq_1_al, std::string& seq_2_al,
	     int d, int* score);

/* sw alignment with affine gap penalty */
int nw_align_ag(std::vector<int>* M, std::vector<int>* I, std::vector<char>* tracebackM,
		std::vector<char>* tracebackI,
		const std::string& seq_1, const std::string& seq_2, 
		std::string& seq_1_al, std::string& seq_2_al,
		int d, int* score, CIGAR_LIST* cigar_list);

void  dpm_init        ( std::vector<int>* F, std::vector<char>* traceback, int L1, int L2, int d);
int   max             ( int f1, int f2, int f3, char * ptr);

int   maxM             ( int f1, int f2, char * ptr);
int   maxI             ( int f1, int f2, int f3, int f4, char * ptr);



