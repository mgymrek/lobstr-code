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
  // Locus properties
  std::string chrom;
  int start;
  int stop;
  std::string repseq;
  std::string repseq_in_ref;
  int period;
  std::string ref_allele;
  float refcopy;
  // Data properties
  int allele1;
  int allele2;
  int coverage;
  float posterior_prob; // posterior prob. of call
  float max_log_lik; // maximum likelihood
  float max_lik_score; // ML/sum of all likelihoods
  float allele1_marginal_lik_score; // marginal likelihood score
  float allele2_marginal_lik_score; // marginal likelihood score
  float allele1_marginal_posterior_prob; // marginal posterior prob
  float allele2_marginal_posterior_prob; // marginal posterior prob
  int conflicting;
  int agreeing;
  int num_stitched;
  std::string readstring;
  std::string allele1_string;
  std::string allele2_string;
  std::map<pair<int,int>,float> likelihood_grid; // log10(P(R|G))
  std::map<int,int> spanning_reads;
  vector<int> alleles_to_include;
  void Reset() {
    chrom = "";
    start = -1;
    stop = -1;
    repseq = "";
    repseq_in_ref = "";
    period = -1;
    ref_allele = "N";
    refcopy = 0;
    allele1 = -10000;
    allele2 = -10000;
    coverage = 0;
    posterior_prob = -1;
    max_log_lik = -10000;
    max_lik_score = -1;
    allele1_marginal_lik_score = -1;
    allele2_marginal_lik_score = -1;
    allele1_marginal_posterior_prob = -1;
    allele2_marginal_posterior_prob = -1;
    conflicting = 0;
    agreeing = 0;
    num_stitched = 0;
    readstring = "";
    allele1_string = "";
    allele2_string = "";
    likelihood_grid.clear();
    spanning_reads.clear();
    alleles_to_include.clear();
  }
};

#endif  // SRC_STRRECORD_H_
