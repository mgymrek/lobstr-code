/*
 * Author: Melissa Gymrek 2012
 */

#ifndef NOISE_MODEL_H_
#define NOISE_MODEL_H_

#include "ReadContainer.h"

using namespace std;

/*
  Class to fit PCR stutter noise
  Can be trained from aligned reads
  or can load a previously trained model
  from a file
 */

class NoiseModel {
 public:
  NoiseModel();
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
  bool HasUniqueMode(list<BamTools::BamAlignment> aligned_reads,
		     float* mode);

  /* Fit logistic model for noise/no noise decision */
  void FitMutationProb();

  /* Fit Poisson model for number of noise steps */
  void FitStepProb();
};

#endif /* NOISE_MODEL_H_ */
