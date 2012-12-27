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

#ifndef SRC_STRRECORD_H_
#define SRC_STRRECORD_H_

using namespace std;

/*
  Struct to keep track of info for a single STR locus
 */
struct STRRecord {
  std::string chrom;
  int start;
  int stop;
  std::string repseq;
  int period;
  int allele1;
  int allele2;
  int coverage;
  float score;
  float allele1_score;
  float allele2_score;
  int conflicting;
  int agreeing;
  int partial_coverage;
  int num_stitched;
  float refcopy;
  int max_partial;
  std::string readstring;
  std::string partialreadstring;
  std::string max_partial_string;
  std::string allele1_string;
  std::string allele2_string;
  std::map<pair<int,int>,float> likelihood_grid;
  std::map<pair<int,int>,float> posterior_grid;
  char ref_allele;
  vector<int> alleles_to_include;
  void Reset() {
    chrom = "";
    start = -1;
    stop = -1;
    repseq = "";
    period = -1;
    allele1 = -10000;
    allele2 = -10000;
    coverage = 0;
    allele1_score = -1;
    allele2_score = -1;
    conflicting = 0;
    agreeing = 0;
    partial_coverage = 0;
    num_stitched = 0;
    refcopy = -1;
    max_partial = 0;
    ref_allele = 'N';
    likelihood_grid.clear();
    posterior_grid.clear();
    alleles_to_include.clear();
  }
};

#endif  // SRC_STRRECORD_H_
