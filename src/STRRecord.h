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
  std::string name;
  // list of samples
  vector<std::string> samples;
  int numcalls;
  // Data properties. each is a vector, indexed by sample number
  vector<int> allele1;
  vector<int> allele2;
  vector<int> coverage;
  vector<float> max_log_lik; // maximum likelihood
  vector<float> max_lik_score; // ML/sum of all likelihoods
  vector<float> allele1_marginal_lik_score; // marginal likelihood score
  vector<float> allele2_marginal_lik_score; // marginal likelihood score
  vector<int> conflicting;
  vector<int> agreeing;
  vector<int> num_stitched;
  vector<std::string> readstring;
  vector<std::string> allele1_string;
  vector<std::string> allele2_string;
  vector<std::map<pair<int,int>,float> > likelihood_grid; // log10(P(R|G))
  vector<std::map<int,int> > spanning_reads;
  vector<int> alleles_to_include;
  vector<float> prob_ref; // P(R|0)/sum all likelihoods
  vector<float> strand_bias;
  vector<float> mean_dist_ends;
  void Reset() {
    chrom = "";
    start = -1;
    stop = -1;
    repseq = "";
    repseq_in_ref = "";
    period = -1;
    ref_allele = "N";
    refcopy = 0;
    samples.clear();
    name = "";
    numcalls = 0;
    allele1.clear();
    allele2.clear();
    coverage.clear();
    prob_ref.clear();
    max_log_lik.clear();
    max_lik_score.clear();
    allele1_marginal_lik_score.clear();
    allele2_marginal_lik_score.clear();
    conflicting.clear();
    agreeing.clear();
    num_stitched.clear();
    readstring.clear();
    allele1_string.clear();
    allele2_string.clear();
    likelihood_grid.clear();
    alleles_to_include.clear();
  }
};

#endif  // SRC_STRRECORD_H_
