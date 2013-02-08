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
#include "src/GenotypeTabWriter.h"
#include "src/runtime_parameters.h"
#include "src/TextFileReader.h"
#include "src/TextFileWriter.h"
#include "src/VCFWriter.h"

using namespace std;
const float SMALL_CONST = 1e-10;
const int PADK = 2;
const float DEFAULT_PRIOR = 0.001;
const float PRIOR_PSEUDOCOUNT = 0.001;

Genotyper::Genotyper(NoiseModel* _noise_model,
                     const vector<string>& _haploid_chroms,
                     map<pair<string,int>, char>* _ref_nucleotides,
                     bool _simple) {
  noise_model = _noise_model;
  haploid_chroms = _haploid_chroms;
  simple = _simple;
  ref_nucleotides = _ref_nucleotides;
  use_priors = false;
}

Genotyper::~Genotyper() {}

void Genotyper::LoadPriors(const std::string& filename) {
  use_priors = true;
  TextFileReader priorfile(filename);
  string line;
  while (priorfile.GetNextLine(&line)) {
    // get locus
    vector<string> items;
    split(line, '\t', items);
    if (items.size() == 0) break;
    if (items.size() != 4) {
      errx(1, "Error, prior file has incorrect number of columns");
    }
    string chrom = items[0];
    int start = atoi(items[1].c_str());
    pair<string, int> locus = pair<string,int>(chrom, start);
    // get allele freqs for each
    map<int, float> allele_to_freq;
    string allele_freq_string = items[3];
    vector<string> alleles;
    split(allele_freq_string, ';', alleles);
    float total = 0;
    for (size_t i = 0; i < alleles.size(); i++) {
      vector<string> astring;
      split(alleles.at(i), ':', astring);
      int allele = atoi(astring[0].c_str());
      float freq = atof(astring[1].c_str());
      allele_to_freq[allele] = freq+PRIOR_PSEUDOCOUNT;
      total += freq;
    }
    // normalize to make sure we have frequencies
    for (map<int, float>::iterator it = allele_to_freq.begin();
         it != allele_to_freq.end(); it++) {
      allele_to_freq[it->first] = it->second/total;
    }
    // insert into allele freqs
    allele_frequencies_per_locus[locus] = allele_to_freq;
  }
  return;
}

float Genotyper::CalcLogLik(int a, int b,
                            const list<AlignedRead>& aligned_reads,
                            int period, int* counta, int* countb ) {
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
    loglik += log10(toadd + SMALL_CONST); // 12/30/31 changed to log10
  }
  return loglik;
}

void Genotyper::FindMLE(const list<AlignedRead>& aligned_reads,
                        const map<int, float>& prior_freqs,
                        bool haploid,
                        int period, int* allele1,
                        int* allele2, float* score, float* maxloglik,
                        float* score_allele1, float* score_allele2,
                        map<pair<int,int>,float>* allelotype_likelihood_grid,
                        map<pair<int,int>,float>* allelotype_posterior_grid,
                        vector<int>* alleles_to_include) {
  int num_nonpartial = 0;
  // Get all possible alleles
  set<int> possible;
  for (list<AlignedRead>::const_iterator it = aligned_reads.begin();
       it != aligned_reads.end(); ++it) {
    if ((*it).partial == 0) {
      possible.insert((*it).diffFromRef);
      num_nonpartial++;
    }
  }
  if (num_nonpartial == 0) {
    *score = -1;
    *score_allele1 = -1;
    *score_allele2 = -1;
    return;
  }

  // Do grid search over (min-k*period, max+k*period)
  map<pair<int,int>,float> allelotype_to_hetfreq;
  int min_allele = *min_element(possible.begin(), possible.end()) - PADK*period;
  int max_allele = *max_element(possible.begin(), possible.end()) + PADK*period;
  // update range based on priors if available
  if (!prior_freqs.empty()) {
    int min_prior_allele = 100000;
    int max_prior_allele = -100000;
    for (map<int,float>::const_iterator it = prior_freqs.begin();
         it != prior_freqs.end(); it++) {
      if (it->first < min_prior_allele) {
        min_prior_allele = it->first;
      }
      if (it->first > max_prior_allele) {
        max_prior_allele = it->first;
      }
    }
    if (min_prior_allele < min_allele) {
      min_allele = min_prior_allele;
    }
    if (max_prior_allele > max_allele) {
      max_allele = max_prior_allele;
    }
  }

  // Get all likelihoods, keep track of max
  float sum_all_posterior = SMALL_CONST; // sum P(R|G)P(G)
  pair<int,int> ml_allelotype(NULL,NULL);
  for (int candidA = min_allele; candidA <= max_allele; candidA++) {
    for (int candidB = candidA; candidB <= max_allele; candidB++) {
      int counta = 0;
      int countb = 0;
      float currScore = CalcLogLik(candidA, candidB,
                                   aligned_reads, period,
                                   &counta, &countb);
      if (debug || plot_info) {
        // Plot info for likelihood grid
        cerr << "[FindMLE]: likelihoods " << aligned_reads.front().chrom << ":" << aligned_reads.front().msStart << " "
             << candidA << " " << candidB << " " << currScore << endl;
      }

      // if haploid and a != b, don't include this as a possible allelotype
      bool include_score = true;
      if ((candidA != candidB) & haploid) {
        include_score = false;
      }
      if (include_score) {
        pair<int, int> allelotype(candidA, candidB);
        allelotype_likelihood_grid->insert(pair<pair<int,int>,float>
                                           (allelotype, currScore));
        float hetfreq = 1.0;
        if (candidA != candidB && !(counta == 0 && countb == 0)) {
          hetfreq = static_cast<float>(counta)/static_cast<float>(counta+countb);
          if (hetfreq > 0.5) hetfreq = 1.0-hetfreq;
        }
        allelotype_to_hetfreq.insert(pair<pair<int,int>,float>
                                     (allelotype, hetfreq));
      }

      // But, include its score in the denominator
      // add P(R|G)P(G) = 10^(log(P(R|G))+log(P(G)))
      float pg = GetPrior(candidA, candidB, prior_freqs);
      sum_all_posterior += pow(10,(currScore+log10(pg)));

      // update which alleles to include in ALT field
      if (include_all_alleles) {
        if (candidA != 0 && find(alleles_to_include->begin(),
                                 alleles_to_include->end(),
                                 candidA) == alleles_to_include->end()
            && !(use_priors && prior_freqs.find(candidA) == prior_freqs.end())) alleles_to_include->push_back(candidA);
        if (candidB != 0 && find(alleles_to_include->begin(),
                                 alleles_to_include->end(),
                                 candidB) == alleles_to_include->end()
            && !(use_priors && prior_freqs.find(candidB) == prior_freqs.end())) alleles_to_include->push_back(candidB);
      }

      // update ML allelotype
      // need to be including score. if including all alleles, must be in alleles to include or 0
      if (include_score && (!include_all_alleles |
                            (include_all_alleles &
                             (find(alleles_to_include->begin(), alleles_to_include->end(), candidA)
                              != alleles_to_include-> end() | candidA == 0) &
                             (find(alleles_to_include->begin(), alleles_to_include->end(), candidB)
                              != alleles_to_include-> end() | candidB == 0)))) {
        if (currScore > *maxloglik) {
          *maxloglik = currScore;
          ml_allelotype.first = candidA;
          ml_allelotype.second = candidB;
        }
      }
    }
  }

  // Get all posteriors, keep track of max
  float bestScore = 0; // max a posteriori
  pair<int,int> best_allelotype(NULL,NULL);
  map<int,float> allele_scores; // marginal posteriors for alleles
  for (map<pair<int,int>,float>::const_iterator it = allelotype_likelihood_grid->begin();
       it != allelotype_likelihood_grid->end(); it++) {
    const int& candidA = it->first.first;
    const int& candidB = it->first.second;
    float pg = GetPrior(candidA, candidB, prior_freqs);
    // update genotype posterior P(R|G)P(G)/sumP(R|G)P(G) = 10^(log(currscore)+log(prior)-log(sum))
    float currScore = pow(10,it->second + log10(pg) - log10(sum_all_posterior));
    // take care of rounding issues, make sure at most 1
    if (currScore > 1) {
      currScore = 1;
    }
    if (debug || plot_info) {
      // Plot info for posterior grid
      cerr << "[FindMLE]: posteriors " << aligned_reads.front().chrom << ":" << aligned_reads.front().msStart << " "
           << candidA << " " << candidB << " " << currScore << endl;
    }
    if (currScore >= bestScore && allelotype_to_hetfreq[it->first] >= min_het_freq) {
      bestScore = currScore;
      best_allelotype = it->first;
    }
    (*allelotype_posterior_grid)[it->first] = currScore;
    // update allele score
    if (allele_scores.find(candidA) == allele_scores.end()) allele_scores[candidA] = 0;
    if (allele_scores.find(candidB) == allele_scores.end()) allele_scores[candidB] = 0;
    allele_scores[candidA] += currScore;
    if (allele_scores[candidA] > 1) allele_scores[candidA] = 1;
    if (candidA != candidB) {
      allele_scores[candidB] += currScore;
      if (allele_scores[candidB] > 1) allele_scores[candidB] = 1;
    }
  }
  if (debug || plot_info) {
    for (map<int,float>::const_iterator it = allele_scores.begin();
         it != allele_scores.end(); it++) {
      cerr << "[FindMLE]: alleleposteriors " << aligned_reads.front().chrom << ":" << aligned_reads.front().msStart << " "
           << it->first << " " << it->second << endl;
    }
  }

  // Return scores
  *allele1 = ml_allelotype.first;
  *allele2 = ml_allelotype.second;
  if (!include_all_alleles && ml_allelotype.first != 0) alleles_to_include->push_back(ml_allelotype.first);
  if (!include_all_alleles &&
      ml_allelotype.first != ml_allelotype.second &&
      ml_allelotype.second != 0) {
    alleles_to_include->push_back(ml_allelotype.second);
  }
  *score = bestScore;
  *score_allele1 = allele_scores[*allele1];
  *score_allele2 = allele_scores[*allele2];
}

void Genotyper::SimpleGenotype(const list<AlignedRead>& aligned_reads,
                               int period,
                               int* allele1, int* allele2,
                               float* score) {
  bool is_haploid = false;
  if ((find(haploid_chroms.begin(), haploid_chroms.end(),
            aligned_reads.front().chrom) != haploid_chroms.end() ||
       find(haploid_chroms.begin(), haploid_chroms.end(), "all")!= haploid_chroms.end())) {
    is_haploid = true;
  }
  if (aligned_reads.size() == 0) {
    *allele1 = 0;
    *allele2 = 0;
    *score = 0;
    return;
  }
  map<float, int> str_to_counts;
  int num_nonpartial = 0;
  for (list<AlignedRead>::const_iterator it = aligned_reads.begin();
       it != aligned_reads.end(); it++) {
    if ((*it).partial == 0 && !(*it).mate) {
      num_nonpartial += 1;
      if (str_to_counts.find((*it).diffFromRef) ==
          str_to_counts.end()) {
        str_to_counts.insert(pair<float, int>((*it).diffFromRef, 1));
      } else {
        str_to_counts.at((*it).diffFromRef) =
          str_to_counts.at((*it).diffFromRef)+1;
      }
    }
  }
  if (num_nonpartial == 0) {
    *allele1 = 0;
    *allele2 = 0;
    *score = 0;
    return;
  }
  int top_str_count = 0;
  float top_str_allele;
  int second_str_count = 0;
  float second_str_allele = 0;
  for (map<float, int>::const_iterator it = str_to_counts.begin();
       it != str_to_counts.end(); it++) {
    float allele = it->first;
    int count = it->second;
    if (count > top_str_count) {
      second_str_count = top_str_count;
      second_str_allele = top_str_allele;
      top_str_count = count;
      top_str_allele = allele;
    } else if (count > second_str_count) {
      second_str_count = count;
      second_str_allele = allele;
    }
  }
  if (is_haploid ||
      static_cast<float>(second_str_count)/
      static_cast<float>(num_nonpartial)
      < min_het_freq) {
    *allele1 = top_str_allele;
    *allele2 = top_str_allele;
    *score = static_cast<float>(top_str_count)/
      static_cast<float>(num_nonpartial);
  } else {
    *allele1 = top_str_allele;
    *allele2 = second_str_allele;
    *score = static_cast<float>(top_str_count + second_str_count) /
      static_cast<float>(num_nonpartial);
  }
}

float Genotyper::GetPrior(int allele1, int allele2,
                          const map<int, float>& prior_freqs) {
  if (prior_freqs.empty()) {
    return 1;
  } else {
    if (prior_freqs.find(allele1) != prior_freqs.end() &&
        prior_freqs.find(allele2) != prior_freqs.end()) {
      return prior_freqs.at(allele1)*prior_freqs.at(allele2);
    } else {
      return DEFAULT_PRIOR;
    }
  }
}

bool Genotyper::ProcessLocus(const std::string chrom,
                             const int str_coord,
                             const std::list<AlignedRead>& aligned_reads,
                             STRRecord* str_record) {
  // Deal with haploid
  bool is_haploid = false;
  if ((find(haploid_chroms.begin(), haploid_chroms.end(), chrom) != haploid_chroms.end() ||
       find(haploid_chroms.begin(), haploid_chroms.end(), "all")!= haploid_chroms.end())) {
    is_haploid = true;
  }

  // Get STR properties
  str_record->coverage = aligned_reads.size();
  if (str_record->coverage == 0) return false;
  str_record->period = aligned_reads.front().period;
  str_record->chrom = chrom;
  str_record->start = aligned_reads.front().msStart;
  str_record->stop = aligned_reads.front().msEnd;
  str_record->repseq = aligned_reads.front().repseq;
  str_record->refcopy = aligned_reads.front().refCopyNum;
  if (ref_nucleotides->find(pair<string,int>(str_record->chrom,
                                             str_record->start)) != ref_nucleotides->end()) {
    str_record->ref_allele = ref_nucleotides->at(pair<string,int>(str_record->chrom,
                                                                  str_record->start));
  }
  if (debug || plot_info) {
    cerr << "[ProcessLocus]: " << str_record->chrom << " " << 
      str_record->start << " " << str_record->repseq << endl;
  }
  // Get alleles and scores
  if (str_record->repseq.empty()) return false;
  if (simple) {
    SimpleGenotype(aligned_reads, str_record->period,
                   &(str_record->allele1), &(str_record->allele2),
                   &(str_record->score));
  } else {
    map<int, float> prior_freqs;
    if (use_priors ) {
      pair<string, int> locus(str_record->chrom, str_record->start);
      if (allele_frequencies_per_locus.find(locus) != allele_frequencies_per_locus.end()) {
        prior_freqs = allele_frequencies_per_locus[locus]; 
      } 
    }
    FindMLE(aligned_reads, prior_freqs, is_haploid, str_record->period,
            &(str_record->allele1), &(str_record->allele2),
            &(str_record->score), &(str_record->max_log_lik),
            &(str_record->allele1_score),
            &(str_record->allele2_score),
            &(str_record->likelihood_grid), 
            &(str_record->posterior_grid),
            &(str_record->alleles_to_include));
    if (debug) {
      cerr << "[ProcessLocus]: mle " << str_record->allele1
           << " " << str_record->allele2 << endl;
    }
  }
  // get read strings
  if (debug) {
    cerr << "[ProcessLocus]: get read strings" << endl;
  }
  map<int, int> all_reads;
  map<int, int> partial_reads;
  float allele = 0;
  bool partial = false;
  int max_partial = -10000;
  for (list<AlignedRead>::const_iterator it2 =
         aligned_reads.begin(); it2 != aligned_reads.end(); it2++) {
    if (it2->mate) continue;
    allele = it2->diffFromRef;
    partial = it2->partial;
    if (it2->stitched && !partial) str_record->num_stitched++;
    if (!partial && (allele == str_record->allele1 ||
                     allele == str_record->allele2)) {
      str_record->agreeing++;
    } else if (!partial) {
      str_record->conflicting++;
    }
    if (partial) {
      if (it2->diffFromRef > max_partial) max_partial = it2->diffFromRef;
      partial_reads[it2->diffFromRef]++;
      str_record->partial_coverage++;
    } else {
      all_reads[it2->diffFromRef]++;
    }
  }
  str_record->coverage = str_record->coverage - str_record->partial_coverage;
  if (debug) {
    cerr << "[ProcessLocus]: coverage " << str_record->coverage << endl;
  }
  // set readstring
  stringstream readstring;
  if (str_record->coverage  == 0) {
    readstring << "NA";
  } else {
    for (map<int, int>::const_iterator vi = all_reads.begin();
         vi != all_reads.end(); vi++) {
      if (vi != all_reads.begin()) {
        readstring << ";";
      }
      readstring << vi->first << "|" << vi->second;
    }
  }
  str_record->readstring = readstring.str();
  // set partial readstring
  stringstream partialreadstring;
  if (str_record->partial_coverage == 0) partialreadstring << "NA";
  for (map<int, int>::const_iterator vi = partial_reads.begin();
       vi != partial_reads.end(); vi++) {
    if (vi != partial_reads.begin()) {
      partialreadstring << ";";
    }
    partialreadstring << vi->first << "|" << vi->second;
  }
  str_record->partialreadstring = partialreadstring.str();
  if (debug) {
    cerr << "[ProcessLocus]: partial coverage "
         << str_record->partial_coverage << endl;
  }
  stringstream max_partial_string;
  if (str_record->partial_coverage != 0) {
    max_partial_string << max_partial;
  } else {max_partial_string << "NA";}
  str_record->max_partial_string = max_partial_string.str();

  // set allele string, deal with low scores
  if (debug) {
    cerr << "[ProcessLocus]: set allele string1 " << endl;
  }
  stringstream allele1_string;
  if (str_record->allele1 == MISSING || str_record->allele2 == MISSING) {
    allele1_string << "NA";
  } else {
    if (str_record->score >= MIN_POSTERIOR || is_haploid) {
      allele1_string << str_record->allele1;
    } else if (str_record->allele1_score >= MIN_MARGINAL) {
      allele1_string << str_record->allele1;
    } else {
      allele1_string << "NA";
    }
  }
  if (debug) {
    cerr << "[ProcessLocus]: set allele string2 " << endl;
  }
  stringstream allele2_string;
  if (str_record->allele1 == MISSING || str_record->allele2 == MISSING) {
    allele2_string << "NA";
  } else if (str_record->score >= MIN_POSTERIOR || is_haploid) {
    allele2_string << str_record->allele2;
  } else if (str_record->allele2_score >= MIN_MARGINAL &&
             !str_record->allele1 == str_record->allele2) {
    allele2_string << str_record->allele2;
  } else {
    allele2_string << "NA";
    if (str_record->allele1 == str_record->allele2) {
      str_record->allele2_score = -1;
    }
  }
  str_record->allele1_string = allele1_string.str();
  str_record->allele2_string = allele2_string.str();
  if (debug) {
    cerr << "[ProcessLocus]: done with process locus " << endl;
  }
  return true;
}

void Genotyper::Genotype(const ReadContainer& read_container,
                         const std::string& output_file,
                         const std::string& vcf_file) {
  GenotypeTabWriter gWriter(output_file);
  VCFWriter vcfWriter(vcf_file);
  STRRecord str_record;
  // collect info on each locus
  for (map<pair<string, int>, list<AlignedRead> >::const_iterator
         it = read_container.aligned_str_map_.begin();
       it != read_container.aligned_str_map_.end(); it++) {
    str_record.Reset();
    if (use_chrom.empty() || (!use_chrom.empty() && use_chrom == it->first.first)) {
      if (ProcessLocus(it->first.first, it->first.second,
                       it->second, &str_record)) {
        gWriter.WriteRecord(str_record);
        vcfWriter.WriteRecord(str_record);
      }
    }
  }
}
