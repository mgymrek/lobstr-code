/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef __FFTW_NUC_VECTORS_H__
#define __FFTW_NUC_VECTORS_H__

#include <fftw3.h>

class FFT_NUC_VECTOR
{
private:
	static fftw_plan plan_256;
	static fftw_plan plan_1024;

	fftw_plan plan;

public:
	fftw_complex *in;
	fftw_complex *out;
	double       *out_magnitude;

	size_t vector_size;

	FFT_NUC_VECTOR() ;
	virtual ~FFT_NUC_VECTOR();
	void resize(size_t minimum_size);
	void allocate_vector(size_t new_size);
	void free_vector();
	void pad_zeros_to_end(size_t start);

	void out_complex_to_magnitude();

	void build_plan(size_t data_size);
	void execute();
	void destroy_plan();

	static void initialize_fftw_plans();
};


#endif
