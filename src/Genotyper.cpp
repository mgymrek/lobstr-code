/*
 * Author: Melissa Gymrek 2012
 */

#include <cmath>
#include <limits>
#include <set>
#include <sstream>

#include "Genotyper.h"
#include "runtime_parameters.h"
#include "TextFileWriter.h"

using namespace std;
const int MIN_PERIOD = 2;
const int MAX_PERIOD = 6;
const float SMALL_CONST = 1e-10;
const float min_supp_freq = 0.5;
const float maxvalue = 1.0/0.0;

static int round(float number) {
  return number < 0.0 ? ceil(number - 0.5) : floor(number + 0.5);
}

Genotyper::Genotyper(NoiseModel* _noise_model,
		     bool _male, bool _simple) {
  noise_model = _noise_model;
  male = _male;
  simple = _simple;
}

Genotyper::~Genotyper(){}

float Genotyper::CalcLogLik(int a, int b,
			    const list<AlignedRead>& aligned_reads,
			    int period, int* counta, int* countb ) {
  // check if in range first
  if (debug) {
    cerr << "in calcloglik " << a << " " << b <<  " " << period << endl;
  }
  float logLik = 0;
  for (list<AlignedRead>::const_iterator
	 it = aligned_reads.begin(); it != aligned_reads.end(); it++) {
    if ((*it).partial == 1) continue;
    int diff = (*it).diffFromRef;
    if (diff == a) {*counta = *counta + 1;}
    if (diff == b) {*countb = *countb + 1;}
    float x = log(SMALL_CONST + noise_model->GetTransitionProb(a, diff, period));
    float y = log(SMALL_CONST + noise_model->GetTransitionProb(b, diff, period));
    float toadd = max(x,y); //log(exp(max(x,y)));//log((exp(y)+exp(x))/2);
    if (debug) {
      cerr << "diff: " << diff << " a: " << a << " b: " << b << " " << x << " " << y << " " << toadd<< endl;
    }
    logLik += toadd;
  }
  if (debug) cerr << logLik << " counta " << *counta << " countb " << *countb << endl;
  return logLik;
}

void Genotyper::FindMLE(const list<AlignedRead>& aligned_reads,
			int period, float* allele1,
			float* allele2, float* score) {
  bool is_haploid = ((aligned_reads.front().chrom == "chrX" ||
		      aligned_reads.front().chrom == "chrY" ) &&
		     male);
  
  float refcopy = aligned_reads.front().refCopyNum;
  if (debug) {
    cerr << "In findMLE " << is_haploid << endl;
  }
  // Get all possible alleles
  set<float> possible;
  for (list<AlignedRead>::const_iterator it = aligned_reads.begin();
       it != aligned_reads.end(); ++it) {
    if ((*it).partial == 0) {
      possible.insert((*it).diffFromRef);
    }
  }
  // if only one possible, homozygous
  if (possible.size() == 1) {
    *allele1 = *(possible.begin());
    *allele2 = *(possible.begin());
    *score = 1;
  }
  // iterate over all possible pairs
  // (only allow homozygous if male and chrX or chrY)
  float maf = 0;
  float perc_supp_reads = 0;
  float currBestScore = -1000;
  float nextBestScore = -1000;
  for (set<float>::const_iterator it1 = possible.begin();
       it1 != possible.end(); it1++) {
    for (set<float>::const_iterator it2 = possible.begin();
	 it2 != possible.end(); it2++) {
      if ((*it2 >= *it1) && !(is_haploid && *it2!=*it1)) {
	int candidA = (int)*it1;
	int candidB = (int)*it2;
	int counta = 0;
	int countb = 0;
	if (debug) {
	  cerr << "calling log lik " << candidA << " " << candidB << endl;
	}
	float currScore = CalcLogLik(candidA, candidB,
				     aligned_reads, period,
				     &counta, &countb);
	// get minor allele freq
	maf = candidA == candidB ? 1 : 
	  (float)min(counta,countb)/(float)(counta+countb);
	perc_supp_reads = candidA == candidB ? (float)counta/(float)aligned_reads.size():
	  (float)(counta + countb)/(float)aligned_reads.size();
	if (debug) {
	  cerr << "perc supp reads " << perc_supp_reads << " maf: " << maf << endl;
	}
	if (currScore >= currBestScore &&
	    maf >= min_het_freq && perc_supp_reads >= min_supp_freq) {
	  nextBestScore = currBestScore;
	  currBestScore = currScore;
	  *allele1 = candidA;
	  *allele2 = candidB;
	  *score = currScore;
	} else if (currScore >= nextBestScore) {
	  nextBestScore = currScore;
	}
      }
    }
  }
  // score is log(e^currBestScore/e^nextBestScore)
  *score = (nextBestScore ==-1000)? maxvalue: currBestScore - nextBestScore;
  if (*score != *score) {*score = 0;} // if score is NaN, set to 0
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
    float allele = aligned_reads.front().diffFromRef;
    *allele1 = allele; *allele2 = allele; *score = 1;
    return;
  }
  map<float, int> str_to_counts;
  for (list<AlignedRead>::const_iterator it = aligned_reads.begin();
       it != aligned_reads.end(); it++) {
    if ((*it).partial==1) continue;
    if (str_to_counts.find((*it).diffFromRef) ==
	str_to_counts.end()) {
      str_to_counts.insert(pair<float,int>((*it).diffFromRef,1));
    } else{
      str_to_counts.at((*it).diffFromRef) =
	str_to_counts.at((*it).diffFromRef)+1;
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
  // write header from bam file
  if (user_defined_arguments.size() > 0) {
    gWriter.Write(user_defined_arguments.substr(4,user_defined_arguments.size()-5));
  }
  // write allelotyper params
  gWriter.Write(user_defined_arguments_allelotyper);
  string chrom;
  int start;
  int stop;
  string repseq;
  int period;
  float allele1 = -10000;
  float allele2 = -10000;
  int allele1_bp;
  int allele2_bp;
  int coverage;
  float score;
  int conflicting;
  int agreeing;
  int partial_coverage;
  float refcopy;
  int max_partial;
  for (map<pair<string, int>, list<AlignedRead> >::const_iterator
	 it = read_container.aligned_str_map_.begin();
       it != read_container.aligned_str_map_.end(); it++) {
    coverage = it->second.size();
    // Get alleles and score
    period = it->second.front().period;
    allele1 = -10000;
    allele2 = -10000;
    if ((it->second.front().chrom == "chrX" ||
	 it->second.front().chrom == "chrY" ) &&
	sex_unknown) {
      // don't genotype sex chroms if gender unknown
      continue;
    }
    if (it->second.front().chrom == "chrY" &&
	!male) {
      // don't genotype chrY if female
      continue;
    }

    if (simple) {
      SimpleGenotype(it->second, period, &allele1, &allele2, &score);
    } else {
      FindMLE(it->second, period, &allele1, &allele2, &score);
    }

    // get read string
    agreeing = 0;
    conflicting = 0;
    partial_coverage = 0;
    max_partial = -1000;
    chrom = "";
    map<int, int> all_reads;
    map<int, int> partial_reads;
    for (list<AlignedRead>::const_iterator it2 = 
	   it->second.begin(); it2 != it->second.end(); it2++) {
      chrom = it2->chrom;
      start = it2->msStart;
      stop = it2->msEnd;
      repseq = it2->repseq;
      refcopy = it2->refCopyNum;
      float allele = it2->diffFromRef;
      if (allele == allele1 || allele == allele2) {
	agreeing++;
	if (allele == allele1) allele1_bp = it2->diffFromRef;
	if (allele == allele2) allele2_bp = it2->diffFromRef;
      } else {
	conflicting++;
      }
      if ((*it2).partial) {
	if (it2->diffFromRef > max_partial) max_partial = it2->diffFromRef;
	partial_reads[it2->diffFromRef]++;
	partial_coverage++;
      } else {
	all_reads[it2->diffFromRef]++;
      }
    }
    stringstream readstring;
    if (coverage - partial_coverage == 0) readstring << "NA";
    for (map<int,int>::const_iterator vi = all_reads.begin();
	 vi != all_reads.end(); vi++) {
      if (vi != all_reads.begin()) {
	readstring << "/";
      }
      readstring << vi->first << ":" << vi->second;
    }
    stringstream partialreadstring;
    if (partial_coverage == 0) partialreadstring << "NA";
    for (map<int,int>::const_iterator vi = partial_reads.begin();
	 vi != partial_reads.end(); vi++) {
      if (vi != partial_reads.begin()) {
	partialreadstring << "/";
      }
      partialreadstring << vi->first << ":" << vi->second;
    }
    stringstream max_partial_string;
    if (partial_coverage != 0) max_partial_string << max_partial;
    else {max_partial_string << "NA";}
    
    stringstream allele1_string;
    if (allele1 == -10000 || allele2 == -10000) {allele1_string << "NA";}
    else {allele1_string << allele1;}
    stringstream allele2_string;
    if (allele1 == -10000 || allele2 == -10000) {allele2_string << "NA";}
    else {allele2_string << allele2;}

    // write to file
    stringstream gLine;
    gLine << chrom << "\t"
	  << start << "\t"
	  << stop << "\t"
	  << repseq << "\t"
	  << period << "\t"
	  << refcopy << "\t"
	  << allele1_string.str() << ","
	  << allele2_string.str() << "\t"
	  << coverage - partial_coverage << "\t"
	  << agreeing << "\t"
	  << conflicting - partial_coverage << "\t"
	  << readstring.str() << "\t"
      	  << score << "\t"
	  << partial_coverage << "\t"
	  << max_partial_string.str() << "\t"
	  << partialreadstring.str();
    
    if (print_reads) {
      gLine << "\t";
      for (list<AlignedRead>::const_iterator readit = it->second.begin();
	   readit != it->second.end(); readit++) {
	if (readit == it->second.begin()) gLine << readit->nucleotides << ":"
						<< readit->diffFromRef << ":" << readit->strand;
	else {
	  gLine << "," << readit->nucleotides<< ":" 
		<< readit->diffFromRef << ":"
		<< readit->strand;
	}
      }
    }
    
    if (!((allele1 == -10000 || allele2 == -10000) && partial_coverage == 0) && (stop > 0 && start > 0 && stop > start)) {
      gWriter.Write(gLine.str());
      if (debug) {
	cerr << gLine.str() << endl;
      }
    }
  }
}
