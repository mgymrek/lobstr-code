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

#ifndef SRC_RUNTIME_PARAMETERS_H_
#define SRC_RUNTIME_PARAMETERS_H_

#include <map>
#include <string>
#include <vector>

// global variables that probably shouldn't be global variables

// stores user arguments for lobSTR alignment
// used to store the bam header for allelotyping step
extern std::string user_defined_arguments;

// store allelotyper user parameters
extern std::string user_defined_arguments_allelotyper;

// keep track of mappings between a kmer and its
// canonical MS so we don't have to recompute it all the time
extern std::map<std::string, std::string> canonicalMSTable;

// enums
enum INPUT_TYPE {
  INPUT_FASTA = 0,
  INPUT_FASTQ,
  INPUT_BAM
};

// flags
extern bool my_verbose;
extern bool sam;
extern bool debug;
extern bool genotype_only;
extern bool fastq;
extern bool bam;
extern bool notab;

// threading
extern size_t threads;

// input files
extern std::string input_files_string;
extern std::string input_files_string_p1;
extern std::string input_files_string_p2;
extern INPUT_TYPE input_type;
extern std::string table_file;
extern std::string genome_file;
extern bool paired;
extern bool gzip;

// output files
extern std::string output_prefix;
extern std::string sam_file;

// detection params
extern bool check_next_best;
extern size_t min_read_length;
extern size_t max_read_length;
extern size_t fft_window_size;
extern size_t fft_window_step;
extern float fft_lobe_threshold;
extern float period_energy_threshold;
extern size_t max_period;
extern size_t min_period;
extern size_t max_period_to_try;
extern size_t min_flank_len;
extern size_t max_flank_len;
extern float closeness;
extern float percent_N_discard;
extern float tukey_alpha;
extern bool use_entropy;
extern float entropy_threshold;
extern int entropy_k;
extern size_t extend_flank;

// alignment params
extern bool adjust;
extern bool debug_adjust;
extern int min_sw_score;
extern int max_mapq;
extern int allowed_mismatches;
extern int max_align;
extern int max_diff_ref;
extern bool unit;
extern int extend;
extern int min_length_to_allow_mismatches;
extern std::string index_prefix;
extern int gap_open;
extern int gap_extend;
extern float fpr;
extern bool partial_debug;

// genotyping params
extern std::string bam_file;
extern std::string command;
extern std::string noise_model;
extern bool rmdup;
extern float min_het_freq;
extern bool male;
extern bool sex_unknown;
extern bool sex_set;
extern int min_coverage;
extern std::string aligned_file;
extern bool non_lobstr_file_detected;
extern bool include_flank;
extern bool print_reads;


// anonymization params (deprecated)
extern std::string str_alignment_file;
extern std::string database_file;
extern bool mask_n;
extern float error_rate;
extern bool include_orig_read_start;

// debug
extern bool profile;
extern bool fftw_debug;
extern bool lobe_debug;
extern bool microsatellite_detection_debug;
extern bool period_detection_debug;
extern bool force_noise_y;
extern int force_noise_y_value;
extern bool multithread_debug;
extern bool gst_debug;
extern bool align_debug;
extern bool anonymizer_debug;
extern bool genotyper_debug;
extern bool why_not_debug;
extern bool entropy_debug;

// min length of str region
extern size_t MIN_STR_LENGTH;

// trimming
extern int QUAL_CUTOFF;
extern int QUALITY_CONSTANT;

#endif  // SRC_RUNTIME_PARAMETERS_H_
