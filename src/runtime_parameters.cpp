/*
 Copyright (C) 2011 Melissa Gymrek <mgymrek@mit.edu>
*/

#include "runtime_parameters.h"

// keep track of common ones
std::map<std::string, std::string> canonicalMSTable;

// flags
bool verbose = false;
bool sam = true;
bool debug = false;
bool genotype_only = false;
bool fastq = false;
bool bam = false;

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
int min_read_length = 45;
int max_read_length = 1024;
int fft_window_size = 24;
int fft_window_step = 12;
float fft_lobe_threshold = 3;
float period_energy_threshold = 500;
int max_period = 6;
int min_period = 2;
int max_period_to_try = 8;
int min_flank_len = 8;
int max_flank_len = 25;
float closeness = 0.3;
float percent_N_discard = 0.05;
float tukey_alpha = 0.5;
bool use_entropy = true;
float entropy_threshold = 0.45;
int entropy_k = 2;
int extend_flank = 6;

// alignment params
bool adjust = true;
int min_sw_score = 10;
int allowed_mismatches = -1;
int max_align = 10;
int max_diff_ref = 30;
bool unit = false;
int extend = 150;
int min_length_to_allow_mismatches = 8;
std::string index_prefix = "";
int gap_open = 1;
int gap_extend = 1;
float fpr = 0.01;


// genotyping params
std::string bam_file = "";
std::string command = "";
std::string noise_model = "";
bool rmdup = true;
float min_het_freq = 0.25;
bool male = true;
bool sex_set = false;
int min_coverage = 2;
std::string aligned_file = "";
bool not_lobstr_file_detected = false;

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
