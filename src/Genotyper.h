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

class Genotyper {
 public:
  Genotyper(NoiseModel* _noise_model,
	    bool _male, bool _simple);
  ~Genotyper();

  /* determine allelotypes and write to file */
  void Genotype(const ReadContainer& read_container,
		const std::string& output_file);
 private:
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
  /* the sample is male */
  bool male;

  /* use the simple genotyper with no noise model */
  bool simple;

  /* store the noise model parameters */
  NoiseModel* noise_model;
};

#endif /* NOISE_MODEL_H_ */
