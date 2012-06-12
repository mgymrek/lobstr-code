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

#include <err.h>
#include <fftw3.h>
#include <pthread.h>

#include <cstdlib>

#include "src/common.h"
#include "src/FFT_nuc_vectors.h"

/*
   FFTW & multi-threading
   ======================

   "execute" is the only thread safe function,
   (see http://www.fftw.org/fftw3_doc/Thread-safety.html#Thread-safety)

   and should be used with the "new-array" variaty
   (see http://www.fftw.org/fftw3_doc/New_002darray-Execute-Functions.html#New_002darray-Execute-Functions)

   We build the static plans once, and later on , all threads will use them.
 */

fftw_plan FFT_NUC_VECTOR::plan_256;
fftw_plan FFT_NUC_VECTOR::plan_1024;

void FFT_NUC_VECTOR::initialize_fftw_plans() {
  fftw_complex *temp_in = reinterpret_cast<fftw_complex*>
    (fftw_malloc(sizeof(fftw_complex) * 1024));
  fftw_complex *temp_out = reinterpret_cast<fftw_complex*>
    (fftw_malloc(sizeof(fftw_complex) * 1024));

  plan_256 = fftw_plan_dft_1d(256, temp_in, temp_out,
                              FFTW_FORWARD, FFTW_ESTIMATE);
  plan_1024 = fftw_plan_dft_1d(1024, temp_in, temp_out,
                               FFTW_FORWARD, FFTW_ESTIMATE);

  fftw_free(temp_in);
  fftw_free(temp_out);
}

FFT_NUC_VECTOR::FFT_NUC_VECTOR()
  : vector_size(0),
    in(NULL), out(NULL), out_magnitude(NULL) {}

FFT_NUC_VECTOR::~FFT_NUC_VECTOR() {
  free_vector();
}

void FFT_NUC_VECTOR::resize(size_t minimum_size) {
  if (vector_size < minimum_size)
    allocate_vector(minimum_size);
}

void FFT_NUC_VECTOR::allocate_vector(size_t new_size) {
  free_vector();
  vector_size = new_size;
  // NOTE: fftw_malloc guarentees WORD-ALIGNMENT in memory.
  // make sure (when replacing with vector<fftw_complex>) to verify alignment.
  // it's required by fftw3 inorder to use SIMD
  in = reinterpret_cast<fftw_complex*>
    (fftw_malloc(sizeof(fftw_complex) * new_size));
  out = reinterpret_cast<fftw_complex*>
    (fftw_malloc(sizeof(fftw_complex) * new_size));
  out_magnitude = reinterpret_cast<double*>
    (malloc(sizeof(double) * new_size)); // NOLINT
}

void FFT_NUC_VECTOR::pad_zeros_to_end(size_t start) {
  for (size_t i = start; i < vector_size; ++i) {
    in[i][0] = 0;
    in[i][1] = 0;
  }
}

void FFT_NUC_VECTOR::free_vector() {
  if (in) {
    fftw_free(in);
    in = NULL;
  }
  if (out) {
    fftw_free(out);
    out = NULL;
  }
  if (out_magnitude) {
    free(out_magnitude);
    out_magnitude = NULL;
  }
}

void FFT_NUC_VECTOR::build_plan(size_t data_size) {
  switch (data_size) {
  case 256:
    plan = plan_256;
    break;

  case 1024:
    plan = plan_1024;
    break;

  default:
    errx(1, "Internal error: no FFTW plan for size = %zu", data_size);
  }
}

void FFT_NUC_VECTOR::execute() {
  fftw_execute_dft(plan, in, out);
}

void FFT_NUC_VECTOR::destroy_plan() {
  // fftw_destroy_plan(plan);
}

void FFT_NUC_VECTOR::out_complex_to_magnitude() {
  for (size_t i = 0; i < vector_size; ++i) {
    out_magnitude[i] = magnitude(out[i]);
  }
}
