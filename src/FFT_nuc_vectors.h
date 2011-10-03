/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>

 This file is part of MicroSatelliteDetector.

 MicroSatelliteDetector is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 MicroSatelliteDetector is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with MicroSatelliteDetector.  If not, see <http://www.gnu.org/licenses/>.
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
