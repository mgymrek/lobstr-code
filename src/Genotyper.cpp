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

#include <algorithm>
#include <cmath>
#include <limits>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>

#include "src/Genotyper.h"
#include "src/runtime_parameters.h"
#include "src/TextFileWriter.h"

using namespace std;
const int MIN_PERIOD = 2;
const int MAX_PERIOD = 6;
const float SMALL_CONST = 1e-10;
const int PADK = 2;

Genotyper::Genotyper(NoiseModel* _noise_model,
                     bool _male, bool _simple) {
  noise_model = _noise_model;
  male = _male;
  simple = _simple;
}

Genotyper::~Genotyper() {}

float Genotyper::CalcLogLik(int a, int b,
                            const list<AlignedRead>& aligned_reads,
                            int period, int* counta, int* countb ) {
  float loglik = 0;
  for (list<AlignedRead>::const_iterator
         it = aligned_reads.begin(); it != aligned_reads.end(); it++) {
    if ((*it).partial == 1 || (*it).mate) continue;
    int diff = (*it).diffFromRef;
    if (diff == a) {*counta = *counta + 1;}
    if (diff == b) {*countb = *countb + 1;}
    float x = noise_model->
      GetTransitionProb(a, diff, period);
    float y = noise_model->
      GetTransitionProb(b, diff, period);
    float toadd = (x+y)/2;
    loglik += log(toadd + SMALL_CONST);
  }
  return loglik;
}

void Genotyper::FindMLE(const list<AlignedRead>& aligned_reads,
                        int period, float* allele1,
                        float* allele2, float* score,
                        float* score_allele1, float* score_allele2) {
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
  map<pair<int,int>,float> allelotype_score_grid;
  map<pair<int,int>,float> allelotype_to_hetfreq;
  int min_allele = *min_element(possible.begin(), possible.end()) - PADK*period;
  int max_allele = *max_element(possible.begin(), possible.end()) + PADK*period;
  if (debug) {
    cerr << "[FindMLE]: grid search from " << min_allele
         << " to " << max_allele << " motif: "
         << aligned_reads.front().repseq << endl;
  }
  float sum_all_likelihoods = 0;
  for (int candidA = min_allele; candidA <= max_allele; candidA++) {
    for (int candidB = candidA; candidB <= max_allele; candidB++) {
      if (debug) {
          cerr << "[FindMLE]: calling log lik " << candidA << " " << candidB
               << " total " << num_nonpartial << endl;
      }
      int counta = 0;
      int countb = 0;
      float currScore = CalcLogLik(candidA, candidB,
                                   aligned_reads, period,
                                   &counta, &countb);
      pair<int, int> allelotype(candidA, candidB);
      allelotype_score_grid.insert(pair<pair<int,int>,float>
                                   (allelotype, currScore));
      float hetfreq = 1.0;
      if (candidA != candidB && !(counta == 0 && countb == 0)) {
        hetfreq = static_cast<float>(counta)/static_cast<float>(counta+countb);
        if (hetfreq > 0.5) hetfreq = 1.0-hetfreq;
      }
      allelotype_to_hetfreq.insert(pair<pair<int,int>,float>
                                   (allelotype, hetfreq));
      sum_all_likelihoods += exp(currScore);
      if (plot_info) {
        cerr << "[FindMLE]: Likelihood " << period << " " << candidA << " "
             << candidB << " " << currScore << " "
             << min_allele << " " << max_allele << endl;
      }
    }
  }

  // Get max posterior
  float bestScore = 0;
  pair<int,int> best_allelotype(NULL,NULL);
  for (map<pair<int,int>,float>::const_iterator it = allelotype_score_grid.begin();
       it != allelotype_score_grid.end(); it++) {
    float currScore = exp(it->second)/sum_all_likelihoods;
    const int& candidA = it->first.first;
    const int& candidB = it->first.second;
    if (plot_info) {
      cerr << "[FindMLE]: Posterior " << period << " " << candidA << " "
           << candidB << " " << currScore
           << " " << min_allele << " " << max_allele << endl;
    }
    if (currScore >= bestScore && allelotype_to_hetfreq[it->first] >= min_het_freq) {
      bestScore = currScore;
      best_allelotype = it->first;
    }
  }

  // Get individual allele scores
  float score1 = 0;
  float score2 = 0;
  for (map<pair<int,int>,float>::const_iterator it = allelotype_score_grid.begin();
       it != allelotype_score_grid.end(); it++) {
    if (it->first.first == best_allelotype.first ||
        it->first.second == best_allelotype.first) {
      score1 += exp(it->second);
    }
    if (it->first.first == best_allelotype.second ||
        it->first.second == best_allelotype.second) {
      score2 += exp(it->second);
    }
  }
  score1 = score1/sum_all_likelihoods;
  score2 = score2/sum_all_likelihoods;
  if (debug) {
    cerr << "[FindMLE]: max posterior " << bestScore << " "
         << "marginal " << best_allelotype.first << " " << score1 << " "
         << "marginal " << best_allelotype.second << " " << score2 << endl;
  }

  // Return scores
  *allele1 = best_allelotype.first;
  *allele2 = best_allelotype.second;
  *score = bestScore;
  *score_allele1 = score1;
  *score_allele2 = score2;
}

void Genotyper::SimpleGenotype(const list<AlignedRead>& aligned_reads,
                               int period,
                               float* allele1, float* allele2,
                               float* score) {
  bool is_haploid = ((aligned_reads.front().chrom == "chrX" ||
                      aligned_reads.front().chrom == "chrY") &&
                     male);
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
    if ((*it).partial == 0) {
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

bool Genotyper::ProcessLocus(const std::string chrom,
                             const int str_coord,
                             const std::list<AlignedRead>& aligned_reads,
                             STRRecord* str_record) {
  // Deal with gender
  if (chrom == "chrY" && !male) {
    // don't genotype chrY if female
    return false;
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
    FindMLE(aligned_reads, str_record->period,
            &(str_record->allele1), &(str_record->allele2),
            &(str_record->score), &(str_record->allele1_score),
            &(str_record->allele2_score));
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
  if (static_cast<float>(str_record->agreeing)/
      static_cast<float>(str_record->coverage) <= min_supp_freq) {
    return false;
  }
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
        readstring << "/";
      }
      readstring << vi->first << ":" << vi->second;
    }
  }
  str_record->readstring = readstring.str();
  if (debug) {
    cerr << readstring.str() << endl;
  }
  // set partial readstring
  stringstream partialreadstring;
  if (str_record->partial_coverage == 0) partialreadstring << "NA";
  for (map<int, int>::const_iterator vi = partial_reads.begin();
       vi != partial_reads.end(); vi++) {
    if (vi != partial_reads.begin()) {
      partialreadstring << "/";
    }
    partialreadstring << vi->first << ":" << vi->second;
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
  stringstream allele1_string;
  if (str_record->allele1 == -10000 || str_record->allele2 == -10000) {
    allele1_string << "NA";
  } else {
    if (str_record->score >= MIN_POSTERIOR) {
      allele1_string << str_record->allele1;
    } else if (str_record->allele1_score >= MIN_MARGINAL) {
      allele1_string << str_record->allele1;
    } else {
      allele1_string << "NA";
    }
  }
  stringstream allele2_string;
  if (str_record->allele1 == -10000 || str_record->allele2 == -10000) {
    allele2_string << "NA";
  } else if (str_record->score >= MIN_POSTERIOR) {
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
  return true;
}

void Genotyper::Genotype(const ReadContainer& read_container,
                         const std::string& output_file) {
  TextFileWriter gWriter(output_file);
  // write header from bam file
  if (user_defined_arguments.size() > 0) {
    gWriter.Write(user_defined_arguments.
                  substr(4, user_defined_arguments.size()-5));
  }
  // write allelotyper params
  gWriter.Write(user_defined_arguments_allelotyper);
  STRRecord str_record;
  // collect info on each locus
  for (map<pair<string, int>, list<AlignedRead> >::const_iterator
         it = read_container.aligned_str_map_.begin();
       it != read_container.aligned_str_map_.end(); it++) {
    str_record.Reset();
    if (ProcessLocus(it->first.first, it->first.second,
                     it->second, &str_record)) {
      // write to file
      stringstream gLine;
      gLine << str_record.chrom << "\t"
            << str_record.start << "\t"
            << str_record.stop << "\t"
            << str_record.repseq << "\t"
            << str_record.period << "\t"
            << str_record.refcopy << "\t"
            << str_record.allele1_string << ","
            << str_record.allele2_string << "\t"
            << str_record.coverage << "\t"
            << str_record.agreeing << "\t"
            << str_record.coverage - str_record.agreeing << "\t"
            << str_record.readstring << "\t"
            << str_record.score << "\t"
            << str_record.allele1_score << "\t"
            << str_record.allele2_score << "\t"
            << str_record.partial_coverage << "\t"
            << str_record.max_partial_string << "\t"
            << str_record.partialreadstring << "\t"
            << str_record.num_stitched;
      if (print_reads) {
        gLine << "\t";
        for (list<AlignedRead>::const_iterator readit = it->second.begin();
             readit != it->second.end(); readit++) {
          if (readit == it->second.begin()) {
            gLine << readit->nucleotides << ":"
                  << readit->diffFromRef << ":"
                  << readit->strand;
          } else {
            gLine << "," << readit->nucleotides<< ":"
                  << readit->diffFromRef << ":"
                  << readit->strand;
          }
        }
      }
      if (!((str_record.allele1 == -10000 || str_record.allele2 == -10000) &&
            str_record.partial_coverage == 0)
          && (str_record.stop > 0 && str_record.start > 0 &&
              str_record.stop > str_record.start)) {
        gWriter.Write(gLine.str());
      }
    }
  }
}
