/*
 * Author: Melissa Gymrek 2012
 */

#include <cmath>
#include <set>
#include <sstream>

#include "Genotyper.h"
#include "runtime_parameters.h"
#include "TextFileWriter.h"

using namespace std;
const int MIN_PERIOD = 2;
const int MAX_PERIOD = 6;
const float SMALL_CONST = 1e-100;

Genotyper::Genotyper(NoiseModel* _noise_model,
		     bool _male, bool _simple) {
  noise_model = _noise_model;
  male = _male;
  _simple = simple;
  if (!_simple)
    PrepareTransitionMatrices();
}

Genotyper::~Genotyper(){}

void Genotyper::PrepareTransitionMatrices() {
  if (debug) {
    cerr << "Preparing transition matrix" << endl;
  }
  for (int period = MIN_PERIOD; period <= MAX_PERIOD; period++) {
    if (debug) {
      cerr << "Period " << period << endl;
    }
    vector<float> tm(MAX_STR_LEN*MAX_STR_LEN);
    transition_matrices.insert(pair<int, vector<float> >
			       (period, tm));
    for (int i = 0; i < MAX_STR_LEN; i++) {
      for (int j = 0; j < MAX_STR_LEN; j++) {
	if (debug) {
	  cerr << "i " << i << " j " << j << endl;
	}
	transition_matrices.at(period).at(i*MAX_STR_LEN+j-1) = 
	  log(SMALL_CONST + noise_model->GetTransitionProb(i, j, period));
      }
    }
  }
}

float Genotyper::CalcLogLik(int a, int b,
			    const list<AlignedRead>& aligned_reads,
			    int period, int* counta, int* countb ) {
  // check if in range first
  if (a > MAX_STR_LEN || b > MAX_STR_LEN) return 0;
  vector<float> joint;
  for (int i = 0; i < MAX_STR_LEN; i++) {
    float x = transition_matrices.at(period-2).at(a*MAX_STR_LEN+i-1);
    float y = transition_matrices.at(period-2).at(b*MAX_STR_LEN+i-1);
    joint.push_back(log((exp(y)+exp(x))/2));
  }
  float logLik = 0;
  for (list<AlignedRead>::const_iterator
	 it = aligned_reads.begin(); it != aligned_reads.end(); it++) {
    float allele = ((*it).diffFromRef/period);
    if ((int)allele == a) {*counta++;}
    if ((int)allele == b) {*countb++;}
    logLik += joint.at((int)allele);
  }
  return logLik;
}

void Genotyper::FindMLE(const list<AlignedRead>& aligned_reads,
			int period, float* allele1,
			float* allele2, float* score) {
  bool is_haploid = ((aligned_reads.front().chrom == "chrX" ||
		      aligned_reads.front().chrom == "chrY" ) &&
		     male);
  if (debug) {
    cerr << "In findMLE " << is_haploid << endl;
  }
  // Get all possible alleles
  set<float> possible;
  for (list<AlignedRead>::const_iterator it = aligned_reads.begin();
       it != aligned_reads.end(); ++it) {
    possible.insert((*it).diffFromRef/period);
  }
  // if only one possible, homozygous
  if (possible.size() == 1) {
    *allele1 = *(possible.begin());
    *allele2 = *(possible.begin());
    *score = 1;
  }
  // iterate over all possible pairs
  // (only allow homozygous if male and chrX or chrY
  float maf = 0;
  float currBestScore = -1000;
  float total_score = 0;
  for (set<float>::const_iterator it1 = possible.begin();
       it1 != possible.end(); it1++) {
    for (set<float>::const_iterator it2 = possible.begin();
	 it2 != possible.end(); it2++) {
      if ((*it2 >= *it1) && !(is_haploid && *it2!=*it1)) {
	float candidA = *it1;
	float candidB = *it2;
	int counta = 0;
	int countb = 0;
	float currScore = CalcLogLik((int)candidA, (int)candidB,
				     aligned_reads, period,
				     &counta, &countb);
	total_score += exp(currScore);
	// get minor allele freq
	maf = counta == countb ? 1 : 
	  (float)max(counta,countb)/(float)(counta+countb);
	if (currScore > currBestScore &&
	    maf >= min_het_freq) {
	  currBestScore = currScore;
	  *allele1 = candidA;
	  *allele2 = candidB;
	  *score = currScore;
	} 
      }
    }
  }
  if (total_score == 0) total_score = 1; // shouldn't happen
  *score = exp(*score)/total_score;
}

void Genotyper::SimpleGenotype(const list<AlignedRead>& aligned_reads,
			       int period,
			       float* allele1, float* allele2,
			       float* score) {
  bool is_haploid = ((aligned_reads.front().chrom == "chrX" ||
		      aligned_reads.front().chrom == "chrY" ) &&
		     male);
  if (aligned_reads.size() == 0) {
    *allele1 = 0; *allele2 = 0; *score = 0;
    return;
  }
  if (aligned_reads.size() == 1) {
    float allele = aligned_reads.front().diffFromRef/period;
    *allele1 = allele; *allele2 = allele; *score = 1;
    return;
  }
  map<float, int> str_to_counts;
  for (list<AlignedRead>::const_iterator it = aligned_reads.begin();
       it != aligned_reads.end(); it++) {
    if (str_to_counts.find((*it).diffFromRef/period) == 
	str_to_counts.end()) {
      str_to_counts.insert(pair<float,int>((*it).diffFromRef/period,1));
    } else{
      str_to_counts.at((*it).diffFromRef/period) = 
	str_to_counts.at((*it).diffFromRef/period)+1;
    }
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
  if (debug) {
    cerr << top_str_count << " " << top_str_allele << " " << second_str_count << " " << second_str_allele << endl;
  }
  if (is_haploid || 
      (float)second_str_count/(float)aligned_reads.size() 
      < min_het_freq) {
    *allele1 = top_str_allele;
    *allele2 = top_str_allele;
    *score = (float)top_str_count/(float)aligned_reads.size();
  } else {
    *allele1 = top_str_allele;
    *allele2 = second_str_allele;
    *score = (float)(top_str_count + second_str_count) / 
      (float)aligned_reads.size();
  }
}

void Genotyper::Genotype(const ReadContainer& read_container,
			 const std::string& output_file) {
  if (debug) {
    cerr << "In main genotype function" << endl;
  }
  TextFileWriter gWriter(output_file);
  string chrom;
  int start;
  int stop;
  string repseq;
  int period;
  float allele1;
  float allele2;
  int coverage;
  float score;
  int conflicting;
  int agreeing;
  float refcopy;
  for (map<pair<string, int>, list<AlignedRead> >::const_iterator
	 it = read_container.aligned_str_map_.begin();
       it != read_container.aligned_str_map_.end(); it++) {
    coverage = it->second.size();
    // Get alleles and score
    period = it->second.front().period;
    if (simple) {
      SimpleGenotype(it->second, period, &allele1, &allele2, &score);
    } else {
      FindMLE(it->second, period, &allele1, &allele2, &score);
    }
    // get read string
    agreeing = 0;
    conflicting = 0;
    map<int, int> all_reads;
    for (list<AlignedRead>::const_iterator it2 = 
	   it->second.begin(); it2 != it->second.end(); it2++) {
       chrom = it2->chrom;
      start = it2->msStart;
      stop = it2->msEnd;
      repseq = it2->repseq;
      refcopy = it2->refCopyNum;
      float allele = it2->diffFromRef/(float)period; 
      if (allele == allele1 || allele == allele2) {
	agreeing++;
      } else {
	conflicting++;
      }
      all_reads[it2->diffFromRef]++;
    }
    stringstream readstring;
    for (map<int,int>::const_iterator vi = all_reads.begin();
	 vi != all_reads.end(); vi++) {
      if (vi != all_reads.begin()) {
	readstring << "/";
      }
      readstring << vi->first << ":" << vi->second;
    }
    // write to file
    stringstream gLine;
    gLine << chrom << "\t"
	  << start << "\t"
	  << stop << "\t"
	  << repseq << "\t"
	  << period << "\t"
	  << refcopy << "\t"
	  << allele1 + refcopy << ","
	  << allele2 + refcopy << "\t"
	  << coverage << "\t"
	  << agreeing << "\t"
	  << conflicting << "\t"
	  << readstring.str() << "\t"
	  << score;
    gWriter.Write(gLine.str());
  }
}
