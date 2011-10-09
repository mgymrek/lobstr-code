/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#include "runtime_parameters.h"

// keep track of common ones
std::map<std::string, std::string> canonicalMSTable;

// flags
bool verbose = false;
bool sam = false;
bool debug = false;
bool genotype_only = false;

// threading
size_t threads = 1;

// input files
std::string input_files_string = "";
INPUT_TYPE input_type = INPUT_FASTA;
std::string table_file = "";
std::string genome_file = "";

// output files
std::string output_prefix = "";
std::string sam_file = "";

// detection params
int min_read_length = 50;
int max_read_length = 1000;
int fft_window_size = 24;
int fft_window_step = 12;
float fft_lobe_threshold = 3;
float period_energy_threshold = 100;
int max_period = 8;
int min_period = 1;
int min_flank_len = 10;
int max_flank_len = 25;
float closeness = 0.5;
float percent_N_discard = 0.2;
float tukey_alpha = 0.5;
bool use_entropy = true;
float entropy_threshold = 0.3;
int entropy_k = 2;
int extend_flank = 5; //8

// alignment params
int allowed_mismatches = 0;
int max_align = 10;
int max_diff_ref = 10;
int extend = 100;
int min_length_to_allow_mismatches = 10;
std::string bwa_ref_prefix = "";

// genotyping params
bool rmdup = false;
std::string mu_string;
std::string pi_string;
float pi_0 = 0.80;
float pi_1 = 0.15;
float pi_2 = 0.05;
float mu_0 = 0.99;
float mu_1 = 0.5;
float mu_2 = 0.01;
bool male = true;
// Note to self: couldn't figure out if it was sexist to have male as default...
int min_coverage = 2;
float percent_reads_for_homozygous = 0.95;
std::string aligned_file = "";

// anonymizing params
std::string str_alignment_file = "";
std::string database_file = "";
bool mask_n = false;
float error_rate = 0.0001;

// debugging
bool profile = false;
bool fftw_debug = false;
bool lobe_debug = false;
bool microsatellite_detection_debug= false ;
bool period_detection_debug = false;
bool force_noise_y = false;
int force_noise_y_value = 2;
bool multithread_debug = false;
bool gst_debug = false;
bool align_debug = false;
bool anonymizer_debug = false;
bool genotyper_debug = false;
bool why_not_debug = false;
bool entropy_debug = false;
