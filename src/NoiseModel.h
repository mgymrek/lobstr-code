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

#include <RInside.h>

#include <list>
#include <string>
#include <vector>

#include "src/ReadContainer.h"

using namespace std;

struct RTrainingData {
  Rcpp::NumericVector mode;
  Rcpp::NumericVector period;
  Rcpp::NumericVector copynum;
  Rcpp::NumericVector mutated;
};

struct TrainingData {
  vector<float> mode;
  vector<int> period;
  vector<float> copynum;
  vector<bool> mutated;
};

/*
  Class to fit PCR stutter noise
  Can be trained from aligned reads
  or can load a previously trained model
  from a file
 */

class NoiseModel {
 public:
  explicit NoiseModel(RInside* _r);
  ~NoiseModel();

  /* Train from a set of aligned reads */
  void Train(ReadContainer* read_container);

  /* Read noise model from file */
  bool ReadNoiseModelFromFile(std::string filename);

  /* Write noise model to file */
  bool WriteNoiseModelToFile(std::string filename);

  /* What is the prob of observing STR=b when true value is STR=a*/
  float GetTransitionProb(int a, int b, int period);

 private:
  /* Check if a set of reads has a unique mode */
  bool HasUniqueMode(const list<AlignedRead>& aligned_reads,
                     float* mode);

  /* Fit logistic model for noise/no noise decision */
  void FitMutationProb();

  /* Fit Poisson model for number of noise steps */
  void FitStepProb();

  /* hold the training data */
  TrainingData training_data;
  RTrainingData r_training_data;

  /* embed an R instance */
  RInside* R;

  /* model data */
  float mutIntercept;
  float mutSlope;
  float poisIntercept;
  float poisSlope;
};

#endif  // SRC_NOISEMODEL_H_
