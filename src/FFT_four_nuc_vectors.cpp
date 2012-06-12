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

#include <string>
#include <vector>

#include "src/common.h"
#include "src/FFT_four_nuc_vectors.h"

using namespace std;

void FFT_FOUR_NUC_VECTORS::allocate_vectors(size_t new_size) {
  A.allocate_vector(new_size);
  C.allocate_vector(new_size);
  G.allocate_vector(new_size);
  T.allocate_vector(new_size);
  summed_matrix.resize(new_size);
}

void FFT_FOUR_NUC_VECTORS::free_vectors() {
  A.free_vector();
  C.free_vector();
  G.free_vector();
  T.free_vector();
}

void FFT_FOUR_NUC_VECTORS::pad_zeros_to_end(size_t start) {
  A.pad_zeros_to_end(start);
  C.pad_zeros_to_end(start);
  G.pad_zeros_to_end(start);
  T.pad_zeros_to_end(start);
}

void FFT_FOUR_NUC_VECTORS::resize(size_t minimum_size) {
  A.resize(minimum_size);
  C.resize(minimum_size);
  G.resize(minimum_size);
  T.resize(minimum_size);
  if (summed_matrix.size() < minimum_size) {
    summed_matrix.resize(minimum_size);
  }
}

void FFT_FOUR_NUC_VECTORS::set_nucleotides(const std::string &nuc) {
  /* populate input data */
  for (size_t i = 0; i < nuc.length(); i++) {
    char ch = nuc[i];
    A.in[i][0] = 0;
    A.in[i][1] = 0;
    C.in[i][0] = 0;
    C.in[i][1] = 0;
    G.in[i][0] = 0;
    G.in[i][1] = 0;
    T.in[i][0] = 0;
    T.in[i][1] = 0;

    switch (ch) {
    case 'a':
    case 'A':
      A.in[i][0] = 1;
    break;

    case 'c':
    case 'C':
      C.in[i][0] = 1;
    break;

    case 'g':
    case 'G':
      G.in[i][0] = 1;
    break;

    case 't':
    case 'T':
      T.in[i][0] = 1;
    break;

    case 'n':
    case 'N':
      continue;
    }
  }
}

void FFT_FOUR_NUC_VECTORS::build_plan(size_t data_size) {
  A.build_plan(data_size);
  C.build_plan(data_size);
  G.build_plan(data_size);
  T.build_plan(data_size);
}

void FFT_FOUR_NUC_VECTORS::execute() {
  A.execute();
  C.execute();
  G.execute();
  T.execute();
}

void FFT_FOUR_NUC_VECTORS::destroy_plan() {
  A.destroy_plan();
  C.destroy_plan();
  G.destroy_plan();
  T.destroy_plan();
}

void FFT_FOUR_NUC_VECTORS::out_complex_to_magnitude() {
  A.out_complex_to_magnitude();
  C.out_complex_to_magnitude();
  G.out_complex_to_magnitude();
  T.out_complex_to_magnitude();
}

void FFT_FOUR_NUC_VECTORS::debug_print_input(size_t stop_at_index) {
  cerr << "Input Data Matrix:" << endl;
  cerr << "xx = [" << endl;
  std::streamsize old_precision = std::cerr.precision();
  std::cerr.precision(4);
  std::ios_base::fmtflags old_flags = std::cerr.flags();
  std::cerr.setf(std::ios::fixed, std::ios::floatfield);
  std::streamsize old_width = std::cerr.width();

  for (size_t i = 0; i < stop_at_index; i++) {
    cerr << A.in[i][0] << "\t"
         << C.in[i][0] << "\t"
         << G.in[i][0] << "\t"
         << T.in[i][0] << "\t"
         << " ;" << endl;
  }
  cerr << "]" << endl;

  std::cerr.precision(old_precision);
  std::cerr.flags(old_flags);
  std::cerr.width(old_width);
}

void FFT_FOUR_NUC_VECTORS::debug_print_output_complex(size_t stop_at_index) {
  cerr << "Output Data Matrix (complex):" << endl;
  cerr << "--- Matlab eqv. fft_xx = fft(xx)" << endl;
  cerr << "fft_xx = [" << endl;
  for (size_t i = 0; i < stop_at_index; i++) {
    cerr << fftw_complex_to_string(A.out[i]) << "\t"
         << fftw_complex_to_string(C.out[i]) << "\t"
         << fftw_complex_to_string(G.out[i]) << "\t"
         << fftw_complex_to_string(T.out[i])
         << " ;" << endl;
  }
  cerr << "]" << endl;
}

void FFT_FOUR_NUC_VECTORS::debug_print_output_magnitude(size_t stop_at_index) {
  cerr << "Output Data Matrix (Magnitude):" << endl;
  cerr << "--- Matlab eqv. mag_fft_x=abs(fft_xx)" << endl;
  cerr << "mag_fft_x = [" << endl;
  for (size_t i = 0; i < stop_at_index; i++) {
    cerr << A.out_magnitude[i] << "\t"
         << C.out_magnitude[i] << "\t"
         << G.out_magnitude[i] << "\t"
         << T.out_magnitude[i]
         << " ;" << endl;
  }
  cerr << "]" << endl;
}

void FFT_FOUR_NUC_VECTORS::create_summed_matrix() {
  for (size_t i = 0; i < summed_matrix.size(); i++) {
    summed_matrix[i] = A.out_magnitude[i] + C.out_magnitude[i] +
      G.out_magnitude[i] + T.out_magnitude[i];
  }
}

void FFT_FOUR_NUC_VECTORS::debug_print_summed_matrix() {
  cerr << "Output Data Matrix (Magnitude):" << endl;
  cerr << "--- matlab eqv. summed_matrix = sum(mag_fft_xx,2)" << endl;
  cerr << "summed_matrix = [" << endl;
  int col = 0;
  for (size_t i = 0; i < summed_matrix.size(); i++) {
    cerr << summed_matrix[i];
    col++;
    if (col == 8) {
      col = 0;
      cerr << "..." << endl;
    } else {
      cerr << "\t";
    }
  }
  cerr << "]" << endl;
}

void FFT_FOUR_NUC_VECTORS::
multiply_nuc_matrix_by_vector(const std::vector<double> &v) {
  if (v.size() > A.vector_size) {
    errx(1, "Internal error (%s:%d): v.size() > A.vector_size",
         __FILE__, __LINE__);
  }
  for (size_t i = 0; i < v.size(); ++i) {
    A.in[i][0] *= v[i];
    C.in[i][0] *= v[i];
    G.in[i][0] *= v[i];
    T.in[i][0] *= v[i];
  }
}

