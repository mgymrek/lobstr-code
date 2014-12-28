/*
Copyright (C) 2011-2014 Melissa Gymrek <mgymrek@mit.edu>

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

RunInfo run_info;
FilterCounter filter_counter;

int MIN_PERIOD=1;
int MAX_PERIOD=6;

// keep track of user defined arguments and git version
std::string user_defined_arguments = "# version=lobSTR_"+std::string(_GIT_VERSION) + ";";
std::string user_defined_arguments_allelotyper = "# version=allelotype_"+std::string(_GIT_VERSION) + ";";
PROGRAM program = LOBSTR;

// flags
bool my_verbose = false;
bool debug = false;
bool fastq = false;
bool bam = false;
bool noweb = false;
bool quiet = false;

// threading
size_t threads = 1;

// input files
std::string input_files_string = "";
std::string input_files_string_p1 = "";
std::string input_files_string_p2 = "";
INPUT_TYPE input_type = INPUT_FASTA;
bool paired = false;
bool gzip = false;

// output files
std::string output_prefix = "";
std::string sam_file = "";

// detection params
size_t min_read_length = 45;
size_t max_read_length = 1024;
size_t fft_window_size = 16;
size_t fft_window_step = 4;
size_t min_flank_len = 8;
size_t max_flank_len = 100;
float percent_N_discard = 0.05;
float entropy_threshold = 0.45;

// alignment params
int min_sw_score = 60;
int max_mapq = 100;
int allowed_mismatches = -1;
int max_align = 10;
int max_diff_ref = 50;
bool unit = false;
int extend = 1000;
int min_length_to_allow_mismatches = 30;
int max_hits_quit_aln = 1000;
bool allow_one_flank_align = true;
std::string index_prefix = "";
int gap_open = 1;
int gap_extend = 1;
float fpr = -1;
int max_mismatch = 1;
std::string read_group_sample = "";
std::string read_group_library = "";
bool allow_multi_mappers = false;

// genotyping params
std::string annotation_files_string = "";
std::string bam_files_string = "";
std::string command = "";
std::string noise_model = "";
std::string haploid_chroms_string = "";
std::string strinfofile = "";
std::string use_chrom = "";
bool rmdup = true;
bool include_gl = false;
bool print_reads = false;
float min_het_freq = 0;
int max_matedist = 100000;
int maximal_end_match_window = 15;
int min_border   = 5;
int min_bp_before_indel = 7;
int min_read_end_match = 5;
bool filter_reads_with_n = false;

// debugging
bool align_debug = false;

size_t MIN_STR_LENGTH = 6;

int QUAL_CUTOFF = 10;
int QUALITY_CONSTANT = 33;
int CHUNKSIZE = 1000;
