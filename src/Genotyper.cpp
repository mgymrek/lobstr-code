/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#include <algorithm>
#include <assert.h>
#include <math.h>
#include <sstream>
#include <boost/math/distributions/binomial.hpp>

#include "common.h"
#include "Genotyper.h"
#include "runtime_parameters.h"
#include "TextFileWriter.h"

using namespace std;
Genotyper::Genotyper() {
  // initialize genotype prior
  genotype_prior.insert(pair<int, float>(0, pi_0));
  genotype_prior.insert(pair<int, float>(1, pi_1));
  genotype_prior.insert(pair<int, float>(2, pi_2));
  
  // initialize read count prior
  read_count_prior.insert(pair<int, float>(0, mu_0));
  read_count_prior.insert(pair<int, float>(1, mu_1));
  read_count_prior.insert(pair<int, float>(2, mu_2));
}

Genotyper::~Genotyper(){};

void Genotyper::AddRead(MSReadRecord* read) {
  pair<string, int> coord (read->chrom, read->msStart);
  // check if this STR is already in the dictionary
  if (aligned_str_map_.find(coord) != aligned_str_map_.end()) {
    aligned_str_map_.at(coord).push_back(*read);
  } else {
    list<MSReadRecord> ms_records_list;
    ms_records_list.push_back(*read);
    aligned_str_map_.insert(pair< pair<string, int>, list<MSReadRecord> >
			    (coord, ms_records_list));
  }
}

void Genotyper::RemovePCRDuplicates(const list<MSReadRecord>& str_records,
				    vector<float>* copy_numbers) {
  if (genotyper_debug) {cout << "removing pcr duplicates..." << endl;}
  // Step 1. group reads into groups of PCR duplicates based
  // on start coordinate and length
  map<ReadPosition, vector<MSReadRecord> > duplicate_groups;
  for (list<MSReadRecord>::const_iterator it = str_records.begin();
       it != str_records.end(); ++it) {
    // get the read position
    ReadPosition read_position;
    read_position.start = (*it).lStart;
    read_position.read_length = (*it).nucleotides.size();
    // insert into the map
    if (duplicate_groups.find(read_position) == duplicate_groups.end()) {
      vector<MSReadRecord> ms_read_records;
      duplicate_groups.insert(pair<ReadPosition, vector<MSReadRecord> >
			      (read_position, ms_read_records));
    }
    duplicate_groups.at(read_position).push_back(*it);
  }
 
  // Step 2. for each group, determine the representative copy number
  // Majority vote if there is one
  // Else vote with highest average quality score
  for (map<ReadPosition, vector<MSReadRecord> >::const_iterator it3 = 
	 duplicate_groups.begin(); it3 != duplicate_groups.end(); ++it3) {
    if (genotyper_debug) {cout << "getting copy for pcr dup group" << endl;}
    // get histogram of genotypes, keep track of qualities 
    // to use as tiebreakers
    map<float, vector<string> > copy_num_to_qualities;
    const vector<MSReadRecord>& records = it3->second;
    for (vector<MSReadRecord>::const_iterator it2 = records.begin();
	 it2 != records.end(); ++it2) {
      if (genotyper_debug) {
	cout << "ref " << (*it2).refCopyNum << endl;
	cout << "diff " << (*it2).diffFromRef << endl;
	cout << "period " << (*it2).ms_repeat_best_period << endl;
      }
      float copy_number = (*it2).refCopyNum +
	(*it2).diffFromRef/(*it2).ms_repeat_best_period;
      if (copy_num_to_qualities.find(copy_number) ==
	  copy_num_to_qualities.end()) {
	vector<string> qualities;
	copy_num_to_qualities.insert(pair<float, vector<string> >
				     (copy_number, qualities));
      }
      copy_num_to_qualities.at(copy_number).push_back((*it2).quality_scores);
    }
    // get majority, use quality as tiebreaker
    float majority_vote_copy_number = -1;
    size_t majority_vote_num_supporting_reads = 0;
    float majority_vote_average_quality = 0;
    for (map<float, vector<string> >::const_iterator it4 =
	   copy_num_to_qualities.begin();
	 it4 != copy_num_to_qualities.end(); ++it4) {
      float average_quality = GetAverageQualityScore(it4->second);
      if (genotyper_debug) {cout << "copy number " << it4->first
				 << " num reads " << it4->second.size()
				 << " average quality " << average_quality << endl;}
      if ( (it4->second.size() > majority_vote_num_supporting_reads) ||
	   ( (it4->second.size() == majority_vote_num_supporting_reads) &&
	     (average_quality > majority_vote_average_quality))) {
	if (genotyper_debug) {cout << "setting copy number" << endl;}
	majority_vote_copy_number = it4->first;
	majority_vote_num_supporting_reads = it4->second.size();
	majority_vote_average_quality = average_quality;
      } 
    }
    // add copy number to the list
    if (majority_vote_copy_number > 0) {
      copy_numbers->push_back(majority_vote_copy_number);
    }
  }
}

pair<float, float> Genotyper::GetGenotype(const vector<float>& copy_numbers,
					  bool autosomal, float ref_copy_number,
					  int* num_conflicting_reads) {
  if (genotyper_debug) {cout << "getting genotype..." << endl;}

  // get read counts map:copynumber->number of supporting reads
  if (genotyper_debug) {cout << "getting read counts..." << endl;}
  map<float, int> copy_number_to_read_counts;
  for (vector<float>::const_iterator it = copy_numbers.begin();
	 it != copy_numbers.end(); ++it) {
    if (copy_number_to_read_counts.find(*it) ==
	copy_number_to_read_counts.end()) {
      copy_number_to_read_counts.insert(pair<float, int> (*it, 1));
    } else {
      copy_number_to_read_counts.at(*it) += 1;
    }
  }

  // get top 2 copy_numbers
  if (genotyper_debug) {cout << "getting top 2 copy numbers" << endl;}
  float top_copy_number = -1;
  size_t top_copy_number_count = 0;
  float second_most_copy_number = -1;
  size_t second_most_copy_number_count = 0;
  for (map<float, int>::const_iterator it = copy_number_to_read_counts.begin();
       it != copy_number_to_read_counts.end(); ++it) {
    if (genotyper_debug) {cout << "copy number " << it->first << " has " << it->second << "reads" << endl;}
    if (it->second > top_copy_number_count) {
      second_most_copy_number = top_copy_number;
      second_most_copy_number_count = top_copy_number_count;
      top_copy_number = it->first;
      top_copy_number_count = it->second;
    } else if (it->second > second_most_copy_number_count) {
      second_most_copy_number = it->first;
      second_most_copy_number_count = it->second;
    }
  }

  if (genotyper_debug) {
    cout << "top genotype " << top_copy_number << "(" 
	 << top_copy_number_count << ")" << endl;
    cout << "second most genotype " << second_most_copy_number << "(" 
	 << second_most_copy_number_count << ")" << endl;
  }

  // if not autosomal, just return majority vote
  if (!autosomal) {
    *num_conflicting_reads = copy_numbers.size() -
      copy_number_to_read_counts.at(top_copy_number);
    return pair<float, float> (top_copy_number, top_copy_number);
  }

  // get posterior probability for each of aa, ab, bb
  if (genotyper_debug) {cout << "get posteriors..." << endl;}
  int genotype = -1; // aa = 0, ab = 1, bb = 2
  float best_posterior_prob = -1;
  for (int i = 0; i < 3; ++i) {
    float prob_i = GetPosteriorProb(copy_numbers, i, ref_copy_number);
    if (genotyper_debug) {
      cout << "genotype " << i 
	   << " posterior " << prob_i << endl;
    }
    if (prob_i > best_posterior_prob) {
      genotype = i;
      best_posterior_prob = prob_i;
    }
  }

  // if aa or ab, return top genotypes
  if (genotyper_debug) {cout << "return genotype... " << genotype << endl;}

  // hack when get ab with only 1 read
  if (second_most_copy_number == -1 && genotype == 1) {
    genotype = 2;
  }

  if (genotype == 0) {
    // assert(top_copy_number == ref_copy_number);
    if (genotyper_debug) {cout << "genotype aa, top copy " 
			       << top_copy_number 
			       << " ref copy " 
			       << ref_copy_number << endl;}
    *num_conflicting_reads = copy_numbers.size() -
      copy_number_to_read_counts.at(ref_copy_number);
    return pair<float, float> (ref_copy_number, ref_copy_number);
  }
  else if (genotype == 1) {
    *num_conflicting_reads = copy_numbers.size() -
      copy_number_to_read_counts.at(top_copy_number) -
      copy_number_to_read_counts.at(second_most_copy_number);
    if (top_copy_number == ref_copy_number) {
      return pair<float, float> (ref_copy_number, second_most_copy_number);
    } else {
      return pair<float, float> (ref_copy_number, top_copy_number);
    }
  }
  else { // if bb
    if (genotyper_debug) {
      cout << "genotype 2 " << "% top copy number reads "
	   << (float(top_copy_number_count)/float(copy_numbers.size())) 
	   << " percent_reads_for_hom " << percent_reads_for_homozygous
	   << endl;
    }
    if (float(top_copy_number_count)/float(copy_numbers.size()) >
	percent_reads_for_homozygous) {
      *num_conflicting_reads = copy_numbers.size() - copy_number_to_read_counts.at(top_copy_number);
      return pair<float, float> (top_copy_number, top_copy_number);
    } else {
       *num_conflicting_reads = copy_numbers.size() -
	 top_copy_number_count - second_most_copy_number_count;
      return pair<float, float> (top_copy_number, second_most_copy_number);
    }
  }
}

bool Genotyper::GetGenotypeString(const pair<string,int>& coordinate,
				  const MSReadRecord& repr_read,
				  string* result) {
  if (genotyper_debug) {
    cout << "getting genotype string..." << endl;
  }
  vector<float> copy_numbers;
  // remove PCR duplicates
  if (rmdup) { RemovePCRDuplicates(aligned_str_map_.at(coordinate), &copy_numbers);}
  else {
    for (list<MSReadRecord>::const_iterator it = aligned_str_map_.at(coordinate).begin();
	 it != aligned_str_map_.at(coordinate).end(); ++it) {
      float new_copy_number = (*it).refCopyNum + (*it).diffFromRef/(*it).ms_repeat_best_period;
      if (new_copy_number > 0) {
	copy_numbers.push_back(new_copy_number);
      }
    }
  }
  if (genotyper_debug) {cout << "got copy numbers coverage..."
			     << copy_numbers.size()
			     << " min coverage " 
			     << min_coverage << endl;}
  // check if above required coverage
  if (copy_numbers.size() < min_coverage) {return false;}
  if (genotyper_debug) { cout << "reporting genotype..." << endl;}
  // get genotype
  bool autosomal = true;
  if (coordinate.first == "chrY" ||
      coordinate.first == "chrX" &&
      male)
    autosomal = false;
  int num_conflicting_reads;
  pair<float, float> genotype = GetGenotype(copy_numbers, autosomal, repr_read.refCopyNum, &num_conflicting_reads);
  sort(copy_numbers.begin(), copy_numbers.end());
  // write the output string
  stringstream ss;
  ss << coordinate.first << "\t"
     << coordinate.second << "\t"
     << repr_read.msEnd << "\t"
     << repr_read.msRepeat << "\t"
     << ((repr_read.name.empty()) ? "-" : repr_read.name) << "\t"
     << repr_read.ms_repeat_best_period << "\t"
     << repr_read.refCopyNum << "\t"
     << min(genotype.first, genotype.second) << ","
     << max(genotype.first, genotype.second) << "\t"
     << copy_numbers.size() << "\t"
     << (copy_numbers.size() - num_conflicting_reads) << "\t"
     << num_conflicting_reads << "\t";
  bool first = true;
  for (vector<float>::const_iterator it = copy_numbers.begin(); it != copy_numbers.end(); ++it) {
    if (first) {ss << *it;}
    else { ss << "/" << *it;}
    first = false;
  }
  //  if (!repr_read.name.empty()) {
  //    ss << "\t" << repr_read.name;
  //}
  *result = ss.str();
  return true;
}

void Genotyper::WriteOutput(const string& filename) {
  if (genotyper_debug) {cout << "Writing genotype output..." << endl;}
  string result;
  TextFileWriter gWriter(filename);
  for (map<pair<string, int>, list<MSReadRecord> >::const_iterator
	 it = aligned_str_map_.begin();
       it != aligned_str_map_.end(); ++it) {
    if (GetGenotypeString(it->first, it->second.front(), &result)) {
      if (genotyper_debug) {cout << "Writing result " << result << endl;}
      gWriter.Write(result);
    }
  }
}

float Genotyper::GetPosteriorProb(const vector<float>& copy_numbers,
				  int genotype, float ref_copy_number) {
  // P(G|reads) = P(reads | genotype) P(genotype) / sum P(g)P(reads|g)
  // but denominator the same for all of them, so only get numerator
  // so technically this isn't the posterior as the function name says
  float p_genotype = genotype_prior.at(genotype);
  // get number of reference reads
  int num_ref_reads = 0;
  for (vector<float>::const_iterator it = copy_numbers.begin();
       it != copy_numbers.end(); ++it) {
    if (*it == ref_copy_number) num_ref_reads++;
  }
  
  // get conditional prob
  // p(reads|genotype) = binom(n,k,p)
  double numTrials = copy_numbers.size();
  double probTrial = read_count_prior.at(genotype);
  boost::math::binomial_distribution<> binom(numTrials, probTrial);
  float p_reads_given_genotype = pdf(binom, num_ref_reads);
  if (genotyper_debug) {
    cout << "genotype " << genotype
	 << " prob_genotype " << p_reads_given_genotype
	 << " " << numTrials << " "
	 << probTrial << " " << num_ref_reads
	 << endl;
    cout << "reads given genotype "
	 << p_reads_given_genotype
	 << " genotype " << p_genotype << endl;
  }
  float posterior = p_reads_given_genotype * p_genotype;
  return posterior;
}

void Genotyper::ResetPi(float pi0, float pi1, float pi2) {
  genotype_prior.clear();
  genotype_prior.insert(pair<int, float>(0, pi0));
  genotype_prior.insert(pair<int, float>(1, pi1));
  genotype_prior.insert(pair<int, float>(2, pi2));
}

void Genotyper::ResetMu(float mu0, float mu1, float mu2) {
  read_count_prior.clear();
  read_count_prior.insert(pair<int, float>(0, mu0));
  read_count_prior.insert(pair<int, float>(1, mu1));
  read_count_prior.insert(pair<int, float>(2, mu2));
}

