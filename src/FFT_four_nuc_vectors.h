/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef SRC_FFT_FOUR_NUC_VECTORS_H__
#define SRC_FFT_FOUR_NUC_VECTORS_H__

#include <string>
#include <vector>

#include "src/FFT_nuc_vectors.h"

class FFT_FOUR_NUC_VECTORS {
 public:
  FFT_NUC_VECTOR A, C, G, T;
  std::vector<double> summed_matrix;

  void allocate_vectors(size_t new_size);
  void free_vectors();
  void pad_zeros_to_end(size_t start);
  void resize(size_t minimum_size);

  void set_nucleotides(const std::string &nuc);
  void multiply_nuc_matrix_by_vector(const std::vector<double> &v);

  void build_plan(size_t data_size);
  void execute();
  void destroy_plan();
  void out_complex_to_magnitude();
  void create_summed_matrix();

  void debug_print_input(size_t stop_at_index);
  void debug_print_output_complex(size_t stop_at_index);
  void debug_print_output_magnitude(size_t stop_at_index);
  void debug_print_summed_matrix();
};

#endif  // SRC_FFT_FOUR_NUC_VECTORS_H__

