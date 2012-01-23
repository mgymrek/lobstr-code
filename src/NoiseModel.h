/*
 * Author: Melissa Gymrek 2012
 */

#ifndef NOISE_MODEL_H_
#define NOISE_MODEL_H_

#include "ReadContainer.h"
#include <RInside.h>

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
  NoiseModel(RInside* _r);
  ~NoiseModel();

  /* Train from a set of aligned reads */
  void Train(ReadContainer* read_container);
  
  /* Read noise model from file */
  bool ReadNoiseModelFromFile(std::string filename);

  /* Write noise model to file */
  bool WriteNoiseModelToFile(std::string filename);

  /* What is the prob of observing STR=b when true value is STR=a*/
  float GetTransitionProb(float a, float b, int period);
     
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

#endif /* NOISE_MODEL_H_ */
