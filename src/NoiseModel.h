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

#ifndef SRC_NOISEMODEL_H_
#define SRC_NOISEMODEL_H_

#include "src/linear.h"

#include <algorithm>
#include <list>
#include <string>
#include <vector>

#include "src/ReadContainer.h"
#include "src/runtime_parameters.h"

using namespace std;

/* Metadata for each STR locus */
struct STRINFO {
  std::string chrom;
  int start;
  int end;
  float score;
  float gc;
  float entropy;
};

/*
  Class to fit PCR stutter noise
  Can be trained from aligned reads
  or can load a previously trained model
  from a file
 */

class NoiseModel {
 public:
  NoiseModel(const std::string& strinfofile,
             const std::vector<std::string>& haploid_chroms_);
  ~NoiseModel();
  
  /* Read STR info from file */
  void ReadSTRInfo(const std::string& filename);

  /* Train from a set of aligned reads */
  void Train(ReadContainer* read_container);

  /* Read noise model from file */
  bool ReadNoiseModelFromFile(const std::string& filename);

  /* What is the prob of observing STR=b when true value is STR=a*/
  float GetTransitionProb(int a, int b, int period, int length,
                          float gc, float score);

  /* STR info data */
  std::map<std::pair<std::string, int>, STRINFO> strInfo;

  std::vector<std::string> haploid_chroms;
 private:
  /* Check if a set of reads has a unique mode */
  bool HasUniqueMode(const std::list<AlignedRead>& aligned_reads,
                     float* mode);

  /* Fit logistic model for noise/no noise decision */
  void FitMutationProb(const std::vector<AlignedRead>& reads_for_training);

  /* Read regression problem from filfe */
  void read_problem(const char *filename);

  /* Fit Poisson model for number of noise steps */
  void FitStepProb(const std::map<int, std::map <int,int> >& step_size_by_period);

  /* model data */
  std::string stutter_problem_filename;
  std::string stutter_model_filename;
  std::string stepsize_model_filename;
  // stutter probability problem
  struct problem stutter_prob;
  // logistic regression of prob. of stutter
  model* stutter_prob_model;
  // probability of stutter to increase allele length
  float p_incr;
  // mean step size mod period
  std::vector<float> average_step_size;
  // PDF of step size for each period
  std::vector<std::vector<float> > pdf_model;
};

#endif  // SRC_NOISEMODEL_H_
