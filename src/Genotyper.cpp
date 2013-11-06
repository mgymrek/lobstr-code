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
#include "src/ZippedTextFileReader.h"

using namespace std;
const float SMALL_CONST = 1e-10;
const int PADK = 2;
const float DEFAULT_PRIOR = 0.001;
const float PRIOR_PSEUDOCOUNT = 0.001;

Genotyper::Genotyper(NoiseModel* _noise_model,
                     const vector<string>& _haploid_chroms,
                     map<pair<string,int>, string>* _ref_nucleotides,
                     map<pair<string,int>, string>* _ref_repseq,
		     const string& vcf_file,
		     const vector<string>& _samples,
		     const map<string,string>& _rg_id_to_sample) {
  noise_model = _noise_model;
  haploid_chroms = _haploid_chroms;
  ref_nucleotides = _ref_nucleotides;
  ref_repseq = _ref_repseq;
  samples = _samples;
  rg_id_to_sample = _rg_id_to_sample;
  vcfWriter = new VCFWriter(vcf_file, samples);
}

Genotyper::~Genotyper() {}

void Genotyper::LoadAnnotations(const vector<std::string> annot_files) {
  string line;
  for (size_t i = 0; i < annot_files.size(); i++) {
    const string& vcf_file = annot_files.at(i);
    if (my_verbose) {
      PrintMessageDieOnError("Loading annotations from file " + vcf_file, PROGRESS);
    }
    // First check that the VCF is zipped and indexed
    if (!fexists(vcf_file.c_str())) {
      PrintMessageDieOnError("File " + vcf_file + " does not exist", ERROR);
    }
    if (vcf_file.find(".vcf.gz") == std::string::npos) {
      PrintMessageDieOnError("File " + vcf_file + " is not a zipped VCF file", ERROR);
    }
    const string index_file = vcf_file + ".tbi";
    if (!fexists(index_file.c_str())) {
      PrintMessageDieOnError("VCF file " + vcf_file + " is not indexed", ERROR);
    }
    ZippedTextFileReader vcfReader(annot_files.at(i));
    while (vcfReader.GetNextLine(&line)) {
      if (!line.empty()) {
	if (line[0] != '#') {
	  STRAnnotation annot;
	  // Get locus info
	  vector<string> items;
	  split(line, '\t', items);
	  annot.chrom = items[0];
	  annot.msStart = atoi(items[1].c_str());
	  annot.name = items[2];
	  int ref_len = static_cast<int>(items[3].size());
	  // Get alt alleles
	  string alt_alleles_string = items[4];
	  vector<string> alt_alleles;
	  split(alt_alleles_string, ',', alt_alleles);
	  annot.alleles.resize(alt_alleles.size());
	  for (size_t j = 0; j < alt_alleles.size(); j++) {
	    int diff = static_cast<int>(alt_alleles.at(j).size())-ref_len;
	    annot.alleles.at(j) = diff;
	  }
	  if (my_verbose) {
	    PrintMessageDieOnError("Loading annotation: " + annot.name, PROGRESS);
	  }
	  annotations[pair<string,int>(annot.chrom, annot.msStart)] = annot;
	}
      }
    }
  }
}

bool Genotyper::GetAlleles(const list<AlignedRead>& aligned_reads,
			   vector<int>* alleles) {
  if (aligned_reads.size() == 0) return false;
  alleles->clear();
  alleles->push_back(0);
  for (list<AlignedRead>::const_iterator it = aligned_reads.begin();
       it != aligned_reads.end(); ++it) {
    if (it->mate) continue;
    if (it->diffFromRef != 0 && std::find(alleles->begin(), alleles->end(), it->diffFromRef) == alleles->end()) {
      alleles->push_back(it->diffFromRef);
    }
  }
  return true;
}

void Genotyper::CleanAllelesList(int reflen, vector<int>* alleles) {
  vector<int> alleles_to_keep(0);
  for (vector<int>::const_iterator it = alleles->begin();
       it != alleles->end(); it++) {
    int total_len = reflen + (*it);
    if (total_len > 0) {
      alleles_to_keep.push_back(*it);
    } else {
      if (debug) {
	PrintMessageDieOnError("Attempt to load invalid allele size", WARNING);
      }
    }
  }
  alleles->resize(alleles_to_keep.size());
  *alleles = alleles_to_keep;
}

// TODO don't read same read twice
bool Genotyper::GetReadsPerSample(const list<AlignedRead>& aligned_reads,
				  const vector<string>& samples,
				  const map<string,string>& rg_id_to_sample,
				  vector<list<AlignedRead> >* sample_reads) {
  sample_reads->clear();
  // Get map of sample->index
  map<string,size_t> sample_to_index;
  for (size_t i = 0; i < samples.size(); i++) {
    sample_to_index[samples[i]] = i;
    sample_reads->push_back(list<AlignedRead>(0));
  }
  // Go through each read and add to appropriate list
  for (list<AlignedRead>::const_iterator it = aligned_reads.begin();
       it != aligned_reads.end(); it++) {
    if (rg_id_to_sample.find(it->read_group) !=
	rg_id_to_sample.end()) {
      string sample = rg_id_to_sample.at(it->read_group);
      string readid = it->ID;
      int i = sample_to_index[sample];
      sample_reads->at(i).push_back(*it);
    } else {
      PrintMessageDieOnError("Could not find sample for read group " + it->read_group, ERROR);
    }
  }
  return true;
}

float Genotyper::CalcLogLik(int a, int b,
                            const list<AlignedRead>& aligned_reads,
                            int period, int* counta, int* countb ) {
  *counta = 0;
  *countb = 0;
  float loglik = 0;
  for (list<AlignedRead>::const_iterator
         it = aligned_reads.begin(); it != aligned_reads.end(); it++) {
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
			map<int,int> spanning_reads,
                        bool haploid, STRRecord* str_record) {
  // Get all likelihoods, keep track of max
  float sum_all_likelihoods = SMALL_CONST; // sum P(R|G)
  // Keep track of numerators for marginal likelihood
  map<int,float> marginal_lik_score_numerator; // sum P(R|G) by allele

  // Things we will add to str record vectors
  map<pair<int,int>, float> likelihood_grid;
  int allele1 = MISSING;
  int allele2 = MISSING;
  float ref_log_lik = -1000000;
  float prob_ref = 0;
  float max_log_lik = -1000000;
  float max_lik_score = 0;
  float allele1_marginal_lik_score = 0;
  float allele2_marginal_lik_score = 0;
  int coverage = aligned_reads.size();
  int conflicting = 0;
  int agreeing = 0;

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
      likelihood_grid.insert
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
      if (debug) {
	stringstream msg;
	msg << "[Genotyper.cpp]: " << currScore << " " << max_log_lik << " "
	    << allelotype.first << "," << allelotype.second;
	PrintMessageDieOnError(msg.str(), DEBUG);
      }
      if (include_score && currScore > max_log_lik &&
          hetfreq >= min_het_freq) {
        max_log_lik = currScore;
        allele1 = allelotype.first;
        allele2 = allelotype.second;
      }
      if (allele1 == 0 && allele2 == 0) {
	ref_log_lik = currScore;
      }
    }
  }

  // Get scores
  max_lik_score =
    pow(10,max_log_lik)/sum_all_likelihoods;
  allele1_marginal_lik_score =
    marginal_lik_score_numerator[allele1]/
    sum_all_likelihoods;
  allele2_marginal_lik_score =
    marginal_lik_score_numerator[allele2]/
    sum_all_likelihoods;
  prob_ref = ref_log_lik - log10(sum_all_likelihoods);

  // Get agreeing/conflicting
  agreeing = spanning_reads[allele1];
  if (allele1 != allele2) {
    agreeing += spanning_reads[allele2];
  }
  conflicting = coverage - agreeing;

  // Add things to STRRecord
  if (coverage != 0) {
    str_record->numcalls++;
  }
  str_record->allele1.push_back(allele1);
  str_record->allele2.push_back(allele2);
  str_record->coverage.push_back(coverage);
  str_record->prob_ref.push_back(prob_ref);
  str_record->max_log_lik.push_back(max_log_lik);
  str_record->max_lik_score.push_back(max_lik_score);
  str_record->allele1_marginal_lik_score.push_back(allele1_marginal_lik_score);
  str_record->allele2_marginal_lik_score.push_back(allele2_marginal_lik_score);
  str_record->conflicting.push_back(conflicting);
  str_record->agreeing.push_back(agreeing);
  str_record->likelihood_grid.push_back(likelihood_grid);
  return;
}

bool Genotyper::ProcessLocus(const std::list<AlignedRead>& aligned_reads,
                             STRRecord* str_record, bool is_haploid) {
  int num_stitched = 0;
  map<int,int> spanning_reads;
  for (list<AlignedRead>::const_iterator it = aligned_reads.begin();
       it != aligned_reads.end(); it++) {
    const AlignedRead& ar = *it;
    if (ar.stitched) num_stitched++;
    spanning_reads[ar.diffFromRef]++;
  }
  str_record->num_stitched.push_back(num_stitched);

  // Get allelotype call and scores
  FindMLE(aligned_reads, spanning_reads, is_haploid, str_record);
  if (debug) {
    stringstream msg;
    msg << "[Genotyper.cpp]: " << str_record->allele1.back() << "," << str_record->allele2.back();
    PrintMessageDieOnError(msg.str(), DEBUG);
  }

  // Get read strings
  stringstream readstring;
  if (aligned_reads.size()  == 0) {
    readstring << "NA";
  } else {
    for (map<int, int>::const_iterator vi = spanning_reads.begin();
         vi != spanning_reads.end(); vi++) {
      if (vi != spanning_reads.begin()) {
        readstring << ";";
      }
      readstring << vi->first << "|" << vi->second;
    }
  }
  str_record->readstring.push_back(readstring.str());
  // set allele string, deal with low scores
  stringstream allele1_string;
  if (str_record->allele1.back() == MISSING || str_record->allele2.back() == MISSING) {
    allele1_string << "NA";
  } else {
    allele1_string << str_record->allele1.back();
  }
  stringstream allele2_string;
  if (str_record->allele1.back() == MISSING || str_record->allele2.back() == MISSING) {
    allele2_string << "NA";
  } else {
    allele2_string << str_record->allele2.back();
  }
  str_record->allele1_string.push_back(allele1_string.str());
  str_record->allele2_string.push_back(allele2_string.str());
  return true;
}
  
void Genotyper::Genotype(const list<AlignedRead>& read_list) {
  STRRecord str_record;
  str_record.Reset();
  // set samples
  str_record.samples = samples;
  // Pull out the chrom and start_coord
  string chrom = read_list.front().chrom;

  // Deal with haploid
  bool is_haploid = false;
  if ((find(haploid_chroms.begin(), haploid_chroms.end(), chrom) != haploid_chroms.end() ||
       find(haploid_chroms.begin(), haploid_chroms.end(), "all") != haploid_chroms.end())) {
    is_haploid = true;
  }

  // Get STR properties
  if (read_list.size() == 0) return;
  str_record.period = read_list.front().period;
  str_record.chrom = chrom;
  if (!(use_chrom.empty() ||
	(!use_chrom.empty() && use_chrom == read_list.front().chrom))) {return;}
  str_record.start = read_list.front().msStart;
  str_record.stop = read_list.front().msEnd;
  str_record.repseq = read_list.front().repseq;
  str_record.refcopy = static_cast<float>((read_list.front().msEnd-
					   read_list.front().msStart))/
    static_cast<float>(read_list.front().period);
  if (ref_nucleotides->find
      (pair<string,int>(str_record.chrom, str_record.start))
      != ref_nucleotides->end()) {
    str_record.ref_allele = ref_nucleotides->at
      (pair<string,int>(str_record.chrom, str_record.start));
    str_record.repseq_in_ref = ref_repseq->at
      (pair<string,int>(str_record.chrom, str_record.start));
  } else {
    return;
  }
  if (str_record.repseq.empty()) return;

  // Check if in our list of annotations
  bool is_annotated = false;
  if (annotations.find(pair<string,int>(str_record.chrom, str_record.start)) !=
      annotations.end()) {
    is_annotated = true;
    const STRAnnotation annot = annotations[pair<string,int>(str_record.chrom, str_record.start)];
    if (my_verbose) {
      PrintMessageDieOnError("Processing annotated locus " + annot.name, PROGRESS);
    }
    str_record.name = annot.name;
    str_record.alleles_to_include = annot.alleles;
  } else {
    // Determine allele range
    if (!GetAlleles(read_list, &str_record.alleles_to_include)) {return;}
  }

  // Clean alleles list to remove things with length < 0
  CleanAllelesList(str_record.stop - str_record.start,
		   &str_record.alleles_to_include);

  // Divide reads for each sample
  vector<list<AlignedRead> > sample_reads;
  if (!GetReadsPerSample(read_list, samples, rg_id_to_sample, &sample_reads)) {return;}

  // Process each sample
  for (size_t i = 0; i < samples.size(); i++) {
    if (debug) {
      PrintMessageDieOnError("[Genotyper.cpp]: Processing sample " + samples.at(i), DEBUG);
    }
    ProcessLocus(sample_reads.at(i), &str_record, is_haploid);
  }
  // Reset alleles to include based on what we called
  if (!is_annotated) {
    str_record.alleles_to_include.clear();
    for (size_t i = 0; i < str_record.samples.size(); i++) {
      int allele1 = str_record.allele1.at(i);
      int allele2 = str_record.allele2.at(i);
      // allele 1
      if (allele1 != 0) {
	if (std::find(str_record.alleles_to_include.begin(),
		      str_record.alleles_to_include.end(),
		      allele1) == str_record.alleles_to_include.end()) {
	  str_record.alleles_to_include.push_back(allele1);
	}
      }
      // allele 2
      if (allele2 != 0) {
	if (std::find(str_record.alleles_to_include.begin(),
		      str_record.alleles_to_include.end(),
		      allele2) == str_record.alleles_to_include.end()) {
	  str_record.alleles_to_include.push_back(allele2);
	}
      }
    }
  }
  std::sort(str_record.alleles_to_include.begin(),
	    str_record.alleles_to_include.end());
  str_record.alleles_to_include.insert(str_record.alleles_to_include.begin(), 0);
  if (str_record.numcalls > 0) {
    vcfWriter->WriteRecord(str_record);
  }
}
