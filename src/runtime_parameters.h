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

#include "src/RunInfo.h"
#include "src/FilterCounter.h"

// period range for genotyper
extern int MIN_PERIOD;
extern int MAX_PERIOD;

// global variables that probably shouldn't be global variables
extern RunInfo run_info;

// when allelotyping, tracks why filtered reads were omitted
extern FilterCounter filter_counter;

// stores user arguments for lobSTR alignment
// used to store the bam header for allelotyping step
extern std::string user_defined_arguments;

// store allelotyper user parameters
extern std::string user_defined_arguments_allelotyper;

// enums
enum INPUT_TYPE {
  INPUT_FASTA = 0,
  INPUT_FASTQ,
  INPUT_BAM
};

enum PROGRAM {
  LOBSTR = 0,
  ALLELOTYPE
};
extern PROGRAM program;

// flags
extern bool my_verbose;
extern bool debug;
extern bool fastq;
extern bool bam;
extern bool noweb;
extern bool quiet;

// threading
extern size_t threads;

// input files
extern std::string input_files_string;
extern std::string input_files_string_p1;
extern std::string input_files_string_p2;
extern INPUT_TYPE input_type;
extern bool paired;
extern bool gzip;

// output files
extern std::string output_prefix;
extern std::string sam_file;

// detection params
extern size_t min_read_length;
extern size_t max_read_length;
extern size_t fft_window_size;
extern size_t fft_window_step;
extern size_t min_flank_len;
extern size_t max_flank_len;
extern float percent_N_discard;
extern float entropy_threshold;

// alignment params
extern int min_sw_score;
extern int max_mapq;
extern int allowed_mismatches;
extern int max_align;
extern int max_diff_ref;
extern bool unit;
extern int extend;
extern int min_length_to_allow_mismatches;
extern int max_hits_quit_aln;
extern bool allow_one_flank_align;
extern std::string index_prefix;
extern int gap_open;
extern int gap_extend;
extern float fpr;
extern int max_mismatch;
extern std::string read_group_sample;
extern std::string read_group_library;
extern bool allow_multi_mappers;

// genotyping params
extern std::string annotation_files_string;
extern std::string bam_files_string;
extern std::string command;
extern std::string noise_model;
extern std::string haploid_chroms_string;
extern std::string strinfofile;
extern std::string use_chrom;
extern bool rmdup;
extern bool include_flank;
extern bool include_gl;
extern bool print_reads;
extern float min_het_freq;
extern int min_border;
extern int min_bp_before_indel;
extern int min_read_end_match;
extern int max_matedist;
extern int maximal_end_match_window;

// debug
extern bool align_debug;

// min length of str region
extern size_t MIN_STR_LENGTH;

// trimming
extern int QUAL_CUTOFF;
extern int QUALITY_CONSTANT;

// chunk size for allelotyper
extern int CHUNKSIZE;

// Amazon s3 parameters
extern bool using_s3;
extern std::string s3bucket;
extern std::string s3cmd_configfile;
extern bool s3debug;

#endif  // SRC_RUNTIME_PARAMETERS_H_
