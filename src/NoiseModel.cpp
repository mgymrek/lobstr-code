/*
 * Author: Melissa Gymrek 2012
 */

#include "NoiseModel.h"

using namespace std;

NoiseModel::NoiseModel(){}

void NoiseModel::Train(ReadContainer* read_container) {
  // Set training data TODO

  // Fit model
  FitMutationProb();
  FitStepProb();
}

bool NoiseModel::HasUniqueMode(list<BamTools::BamAlignment>
			       aligned_reads,
			       float* mode) {
  // TODO
  return false;
}

void NoiseModel::FitMutationProb() {
  // TODO
}

void NoiseModel::FitStepProb() {
  // TODO
}

float NoiseModel::GetTransitionProb(float a,
				    float b, int period) {
  // TODO
  return 0.0;
}

bool NoiseModel::ReadNoiseModelFromFile(string filename) {
  // TODO
  return false;
}

bool NoiseModel::WriteNoiseModelToFile(string filename) {
  // TODO
}

NoiseModel::~NoiseModel(){}
