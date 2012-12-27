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

#include <map>
#include <string>
#include "src/runtime_parameters.h"

// keep track of user defined arguments
std::string user_defined_arguments = "# version=lobSTR_v2.0.2;";
std::string user_defined_arguments_allelotyper = "# version=allelotype_v2.0.2;";

// keep track of common ones
std::map<std::string, std::string> canonicalMSTable;

// flags
bool my_verbose = false;
bool sam = true;
bool debug = false;
bool plot_info = false;
bool genotype_only = false;
bool fastq = false;
bool bam = false;
bool notab = false;

// threading
size_t threads = 1;

// input files
std::string input_files_string = "";
std::string input_files_string_p1 = "";
std::string input_files_string_p2 = "";
INPUT_TYPE input_type = INPUT_FASTA;
std::string table_file = "";
std::string genome_file = "";
bool paired = false;
bool gzip = false;

// output files
std::string output_prefix = "";
std::string sam_file = "";

// detection params
bool check_next_best = true;
size_t min_read_length = 45;
size_t max_read_length = 1024;
size_t fft_window_size = 24;
size_t fft_window_step = 12;
float fft_lobe_threshold = 3;
float period_energy_threshold = 500;
size_t max_period = 6;
size_t min_period = 2;
size_t max_period_to_try = 6;
size_t min_flank_len = 8;
size_t max_flank_len = 25;
float closeness = 0.3;
float percent_N_discard = 0.05;
float tukey_alpha = 0.5;
bool use_entropy = true;
float entropy_threshold = 0.45;
int entropy_k = 2;
size_t extend_flank = 6;

// alignment params
bool adjust = true;
bool debug_adjust = false;
int min_sw_score = 60;
int max_mapq = 100;
int allowed_mismatches = -1;
int max_align = 10;
int max_diff_ref = 50;
bool unit = false;
int extend = 1000;
int min_length_to_allow_mismatches = 8;
std::string index_prefix = "";
int gap_open = 1;
int gap_extend = 1;
float fpr = 0.01;
bool partial_debug = false;
std::string read_group = "";

// genotyping params
std::string bam_files_string = "";
std::string command = "";
std::string noise_model = "";
bool rmdup = true;
float min_het_freq = 0.20;
bool male = true;
bool sex_unknown = false;
bool sex_set = false;
int min_coverage = 1;
std::string aligned_file = "";
bool not_lobstr_file_detected = false;
bool include_flank = true;
bool print_reads = false;
int max_matedist = 100000;
float min_supp_freq = 0.25;
float MIN_POSTERIOR = 0.0;
float MIN_MARGINAL = 0.0;
std::string haploid_chroms_string = "";
bool exclude_partial = false;
std::string strinfofile = "";
std::string priorsfile = "";

// vcf params
std::string exclude_positions_file = "";
std::string sample = "";

// anonymizing params
std::string str_alignment_file = "";
std::string database_file = "";
bool mask_n = false;
float error_rate = 0.0001;
bool include_orig_read_start = false;

// debugging
bool profile = false;
bool fftw_debug = false;
bool lobe_debug = false;
bool microsatellite_detection_debug = false;
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

size_t MIN_STR_LENGTH = 6;

int QUAL_CUTOFF = 15;
int QUALITY_CONSTANT = 33;

// Amazon s3 paramters
bool using_s3 = false;
std::string s3bucket = "";
std::string s3cmd_configfile = "";
bool s3debug = false;
