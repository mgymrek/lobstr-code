/*
 * Author: Melissa Gymrek 2012
 */

#ifndef GENOTYPER_H_
#define GENOTYPER_H_

#include <vector>

#include "NoiseModel.h"
#include "ReadContainer.h"

using namespace std;

/*
  Class to determine allelotypes at each locus
 */

const int MAX_STR_LEN = 100;

class Genotyper {
 public:
  Genotyper(const NoiseModel& _noise_model,
	    bool _male, bool _simple);
  ~Genotyper();

  /* determine allelotypes and write to file */
  void Genotype(const ReadContainer& read_container,
		const std::string& output_file);
 private:
  /* Prepare the transition matrix */
  void PrepareTransitionMatrices();

  /* Get log likelihood of an allelotype */
  float CalcLogLik(int a, int b,
		   const list<AlignedRead>& aligned_reads,
		   int period, int* counta, int* countb);

  /* Get most likely allelotype */
  void FindMLE(const list<AlignedRead>& aligned_reads, int period,
	       float* allele1, float* allele2, float* score);

  /* Get allelotype without using noise model */
  void SimpleGenotype(const list<AlignedRead>& aligned_reads,
		      int period,
		      float* allele1, float* allele2, float* score);
  bool male;
  bool simple;
  NoiseModel noise_model;

  // store transition matrices for each period
  map<int, vector<float> > transition_matrices;
};

#endif /* NOISE_MODEL_H_ */
