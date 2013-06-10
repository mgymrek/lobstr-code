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

#include <err.h>
#include <errno.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>

#include "src/common.h"
#include "src/Genotyper.h"
#include "src/runtime_parameters.h"
#include "src/TextFileReader.h"
#include "src/TextFileWriter.h"

using namespace std;
const float SMALL_CONST = 1e-10;
const int PADK = 2;
const float DEFAULT_PRIOR = 0.001;
const float PRIOR_PSEUDOCOUNT = 0.001;

Genotyper::Genotyper(NoiseModel* _noise_model,
                     const vector<string>& _haploid_chroms,
                     map<pair<string,int>, string>* _ref_nucleotides,
                     map<pair<string,int>, string>* _ref_repseq,
		     const string& vcf_file) {
  noise_model = _noise_model;
  haploid_chroms = _haploid_chroms;
  ref_nucleotides = _ref_nucleotides;
  ref_repseq = _ref_repseq;

  vcfWriter = new VCFWriter(vcf_file);
}

Genotyper::~Genotyper() {}

float Genotyper::CalcLogLik(int a, int b,
                            const list<AlignedRead>& aligned_reads,
                            int period, int* counta, int* countb ) {
  *counta = 0;
  *countb = 0;
  float loglik = 0;
  for (list<AlignedRead>::const_iterator
         it = aligned_reads.begin(); it != aligned_reads.end(); it++) {
    if ((*it).partial == 1 || (*it).mate) continue;
    int diff = (*it).diffFromRef;
    int length = (*it).msEnd-(*it).msStart + diff;
    pair<string, int> coord((*it).chrom, (*it).msStart);
    float gc = noise_model->strInfo[coord].gc;
    float score = noise_model->strInfo[coord].score;
    if (diff == a) {*counta = *counta + 1;}
    if (diff == b) {*countb = *countb + 1;}
    float x = noise_model->
      GetTransitionProb(a, diff, period, length, gc, score);
    float y = noise_model->
      GetTransitionProb(b, diff, period, length, gc, score);
    float toadd = (x+y)/2;
    loglik += log10(toadd + SMALL_CONST);
  }
  return loglik;
}

void Genotyper::FindMLE(const list<AlignedRead>& aligned_reads,
                        bool haploid, STRRecord* str_record) {
  // Get all possible alleles and set other info while we're at it
  float allele;
  bool partial = false;
  set<int> possible;
  for (list<AlignedRead>::const_iterator it = aligned_reads.begin();
       it != aligned_reads.end(); ++it) {
    if (it->mate) continue;
    partial = it->partial;
    allele = it->diffFromRef;
    if (!partial) {
      possible.insert((*it).diffFromRef);
      str_record->coverage++;
      str_record->spanning_reads[allele]++;
    } else {
      str_record->partial_coverage++;
      str_record->partial_reads[allele]++;
      if (allele > str_record->max_partial) {
        str_record->max_partial = allele;
      }
    }
    if (it->stitched && !partial) str_record->num_stitched++;
  }
  if (str_record->coverage == 0) {
    return;
  }

  // Do grid search over (min-k*period, max+k*period)
  int min_allele = *min_element(possible.begin(), possible.end()) - PADK*str_record->period;
  int max_allele = *max_element(possible.begin(), possible.end()) + PADK*str_record->period;
  // Update alleles to include
  // Add 0 right away at the front
  str_record->alleles_to_include.insert
    (str_record->alleles_to_include.begin(), 0);
  for (int i = min_allele; i <= max_allele; i++) {
    if (i != 0) {
      if (str_record->refcopy * str_record->period > -1*i) {
	str_record->alleles_to_include.push_back(i);
      } else {
	stringstream msg;
	msg << "Attempted to load invalid allele size at "
	    << str_record->chrom << ":" << str_record->start;
	PrintMessageDieOnError(msg.str(), WARNING);
      }
    }
  }

  // Get all likelihoods, keep track of max
  float sum_all_likelihoods = SMALL_CONST; // sum P(R|G)
  // Keep track of numerators for marginal likelihood
  map<int,float> marginal_lik_score_numerator; // sum P(R|G) by allele
  for (size_t i = 0; i < str_record->alleles_to_include.size(); i++) {
    for (size_t j = i; j < str_record->alleles_to_include.size(); j++) {
      pair<int, int> allelotype(str_record->alleles_to_include.at(i),
                                str_record->alleles_to_include.at(j));
      if (allelotype.first > allelotype.second) {
        allelotype.first = str_record->alleles_to_include.at(j);
        allelotype.second = str_record->alleles_to_include.at(i);
      }
      int counta,countb;
      float currScore = CalcLogLik(allelotype.first, allelotype.second,
                                   aligned_reads, str_record->period,
                                   &counta, &countb);
      // Insert likelihood to grid
      str_record->likelihood_grid.insert
        (pair<pair<int,int>,float>(allelotype, currScore));
      // Insert het freq into grid
      float hetfreq = 1.0;
      if (allelotype.first != allelotype.second && !(counta == 0 && countb == 0)) {
        hetfreq = static_cast<float>(counta)/static_cast<float>(counta+countb);
        if (hetfreq > 0.5) hetfreq = 1.0-hetfreq;
      }
      // if haploid and a != b, don't include as a possibility
      bool include_score = true;
      if ((allelotype.first != allelotype.second) & haploid) {
        include_score = false;
      }
      // But, include its score in the denominator
      // add P(R|G)P(G) = 10^(log(P(R|G))+log(P(G)))
      float likelihood_term = pow(10, currScore);
      sum_all_likelihoods += likelihood_term;
      marginal_lik_score_numerator[allelotype.first] += likelihood_term;
      if (allelotype.first != allelotype.second) {
        marginal_lik_score_numerator[allelotype.second] += likelihood_term;
      }

      // update ML allelotype, must have include_score true
      if (include_score && currScore > str_record->max_log_lik &&
          hetfreq >= min_het_freq) {
        str_record->max_log_lik = currScore;
        str_record->allele1 = allelotype.first;
        str_record->allele2 = allelotype.second;
      }
    }
  }

  // Get scores
  str_record->max_lik_score =
    pow(10,str_record->max_log_lik)/sum_all_likelihoods;
  str_record->allele1_marginal_lik_score =
    marginal_lik_score_numerator[str_record->allele1]/
    sum_all_likelihoods;
  str_record->allele2_marginal_lik_score =
    marginal_lik_score_numerator[str_record->allele2]/
    sum_all_likelihoods;

  // Get agreeing/conflicting
  str_record->agreeing =
    str_record->spanning_reads[str_record->allele1];
  if (str_record->allele1 != str_record->allele2) {
    str_record->agreeing += str_record->spanning_reads[str_record->allele2];
  }
  str_record->conflicting =
    str_record->coverage -
    str_record->agreeing;
}

bool Genotyper::ProcessLocus(const std::list<AlignedRead>& aligned_reads,
                             STRRecord* str_record) {
  // Pull out the chrom and start_coord
  string chrom = aligned_reads.front().chrom;

  // Deal with haploid
  bool is_haploid = false;
  if ((find(haploid_chroms.begin(), haploid_chroms.end(), chrom) != haploid_chroms.end() ||
       find(haploid_chroms.begin(), haploid_chroms.end(), "all") != haploid_chroms.end())) {
    is_haploid = true;
  }

  // Get STR properties
  if (aligned_reads.size() == 0) return false;
  str_record->period = aligned_reads.front().period;
  str_record->chrom = chrom;
  str_record->start = aligned_reads.front().msStart;
  str_record->stop = aligned_reads.front().msEnd;
  str_record->repseq = aligned_reads.front().repseq;
  str_record->refcopy = static_cast<float>((aligned_reads.front().msEnd-aligned_reads.front().msStart))/
    static_cast<float>(aligned_reads.front().period);
  if (ref_nucleotides->find
      (pair<string,int>(str_record->chrom, str_record->start))
      != ref_nucleotides->end()) {
    str_record->ref_allele = ref_nucleotides->at
      (pair<string,int>(str_record->chrom, str_record->start));
    str_record->repseq_in_ref = ref_repseq->at
      (pair<string,int>(str_record->chrom, str_record->start));
  } else {
    return false;
  }
  if (str_record->repseq.empty()) return false;
  // Get allelotype call and scores
  FindMLE(aligned_reads, is_haploid, str_record);
  // Get read strings
  stringstream readstring;
  if (str_record->coverage  == 0) {
    readstring << "NA";
  } else {
    for (map<int, int>::const_iterator vi = str_record->spanning_reads.begin();
         vi != str_record->spanning_reads.end(); vi++) {
      if (vi != str_record->spanning_reads.begin()) {
        readstring << ";";
      }
      readstring << vi->first << "|" << vi->second;
    }
  }
  str_record->readstring = readstring.str();
  // set partial readstring
  stringstream partialreadstring;
  if (str_record->partial_coverage == 0) partialreadstring << "NA";
  for (map<int, int>::const_iterator vi = str_record->partial_reads.begin();
       vi != str_record->partial_reads.end(); vi++) {
    if (vi != str_record->partial_reads.begin()) {
      partialreadstring << ";";
    }
    partialreadstring << vi->first << "|" << vi->second;
  }
  str_record->partialreadstring = partialreadstring.str();
  stringstream max_partial_string;
  if (str_record->partial_coverage != 0) {
    max_partial_string << str_record->max_partial;
  } else {max_partial_string << "NA";}
  str_record->max_partial_string = max_partial_string.str();

  // set allele string, deal with low scores
  stringstream allele1_string;
  if (str_record->allele1 == MISSING || str_record->allele2 == MISSING) {
    allele1_string << "NA";
  } else {
    allele1_string << str_record->allele1;
  }
  stringstream allele2_string;
  if (str_record->allele1 == MISSING || str_record->allele2 == MISSING) {
    allele2_string << "NA";
  } else {
    allele2_string << str_record->allele2;
  }
  str_record->allele1_string = allele1_string.str();
  str_record->allele2_string = allele2_string.str();
  return true;
}
  
void Genotyper::Genotype(const list<AlignedRead>& read_list) {
  STRRecord str_record;
  str_record.Reset();
  if (use_chrom.empty() ||
      (!use_chrom.empty() && use_chrom == read_list.front().chrom)) {
    if (ProcessLocus(read_list, &str_record)) {
      vcfWriter->WriteRecord(str_record);
    }
  }
}
