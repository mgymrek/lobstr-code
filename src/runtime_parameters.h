/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#ifndef RUNTIME_PARAMETERS_H_
#define RUNTIME_PARAMETERS_H_

#include <map>
#include <string>
#include <vector>

// global variables that probably shouldn't be global variables

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
extern bool verbose;
extern bool sam;
extern bool debug;
extern bool genotype_only;

// threading
extern size_t threads;

// input files
extern std::string input_files_string;
extern INPUT_TYPE input_type;
extern std::string table_file;
extern std::string genome_file;

// output files
extern std::string output_prefix;
extern std::string sam_file;

// detection params
extern int min_read_length;
extern int max_read_length;
extern int fft_window_size;
extern int fft_window_step;
extern float fft_lobe_threshold;
extern float period_energy_threshold;
extern int max_period;
extern int min_period;
extern int min_flank_len;
extern int max_flank_len;
extern float closeness;
extern float percent_N_discard;
extern float tukey_alpha;
extern bool use_entropy;
extern float entropy_threshold;
extern int entropy_k;
extern int extend_flank;

// alignment params
extern int allowed_mismatches;
extern int max_align;
extern int max_diff_ref;
extern int extend;
extern int min_length_to_allow_mismatches;

// genotyping params
extern bool rmdup;
extern float homozygous_entropy_threshold;
extern std::string mu_string;
extern std::string pi_string;
extern float pi_0;
extern float pi_1;
extern float pi_2;
extern float mu_0;
extern float mu_1;
extern float mu_2;
extern bool male;
extern int min_coverage;
extern float percent_reads_for_homozygous;
extern std::string aligned_file;

// anonymization params
extern std::string str_alignment_file;
extern std::string database_file;
extern bool mask_n;
extern float error_rate;

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

#endif /* RUNTIME_PARAMETERS_H_ */
