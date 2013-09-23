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
#include <getopt.h>
#include <limits.h>
#include <stdlib.h>
#include <unistd.h>

#include <boost/algorithm/string.hpp>
#include <iostream>
#include <list>
#include <map>
#include <string>
#include <vector>
#include <utility>

#include "src/BWAReadAligner.h"
#include "src/bwtaln.h"
#include "src/bwase.h"
#include "src/common.h"
#include "src/FastaFileReader.h"
#include "src/FastqFileReader.h"
#include "src/FFT_four_nuc_vectors.h"
#include "src/FFT_nuc_vectors.h"
#include "src/HammingWindowGenerator.h"
#include "src/IFileReader.h"
#include "src/MSReadRecord.h"
#include "src/MultithreadData.h"
#include "src/SamFileWriter.h"
#include "src/STRDetector.h"
#include "src/runtime_parameters.h"
#include "src/TukeyWindowGenerator.h"

using namespace std;
const int READPROGRESS = 10000;
// list of input files to process from
vector<string> input_files;
vector<string> input_files1;
vector<string> input_files2;

// keep track of # bases so we can calculate coverage
int bases = 0;

// keep track of reference genome to use for sam format
map<string, int> chrom_sizes;

// Keep track of reference sequences for alignment readjustment
map<int, REFSEQ> ref_sequences;

// alignment references, keep global
void LoadReference(const std::string& repseq);
void DestroyReferences();
map<std::string, BWT> bwt_references;
map<std::string, BNT> bnt_annotations;
gap_opt_t *opts;

void show_help() {
  const char* help =
    "\nlobSTR [OPTIONS] " \
    "    {-f <file1[,file2,...]> | --p1 <file1_1[,file2_1,...]>\n" \
    "    --p2 <file1_2[,file2_1,...]>} --index-prefix <index prefix>\n" \
    "    -o <output prefix>\n" \
    "Note: parameters are uploaded to Amazon S3 by default. This is for\n" \
    "us see how people are using the tool and to help us continue to improve\n" \
    "lobSTR. To turn this function off, specify --noweb.\n\n" \
    "Parameter descriptions:\n " \
    "-f,--files     file or comma-separated list of files\n" \
    "               containing reads in fasta, fastq, or bam format\n" \
    "               (default: fasta)\n" \
    "--p1           file or comma-separated list of files containing\n" \
    "               the first end of paired end reads in fasta or fastq\n" \
    "               (default: fasta)\n" \
    "--p2           file or comma-separated list of files containing\n" \
    "               the second end of paired end reads in fasta or fastq\n" \
    "               (default: fasta)\n" \
    "-o,--out       prefix for output files. will output:\n" \
    "                  <prefix>.aligned.bam: bam file of alignments\n" \
    "                  <prefix>.aligned.stats: give statistics about alignments\n" \
    "--index-prefix prefix for lobSTR's bwa reference (must run lobstr_index.py\n" \
    "               to create index. If the index is downloaded\n" \
    "               to PATH_TO_INDEX, this argument is\n" \
    "               PATH_TO_INDEX/lobSTR_)\n" \
    "\n\nOptions:\n" \
    "-h,--help      display this help screen\n" \
    "-v,--verbose   print out useful progress messages\n" \
    "--version      print out lobSTR program version\n" \
    "-q,--fastq     reads are in fastq format (default: fasta)\n" \
    "--bam          reads are in bam format (default: fasta)\n" \
    "--gzip         The input files are gzipped\n" \
    "               (only works for fasta or fastq input)\n" \
    "--bampair      reads are in bam format and are paired-end\n" \
    "               NOTE: bam file MUST be sorted by read name\n" \
    "               (samtools sort -n <file.bam> <prefix>)\n" \
    "--bwaq         Trim read ends based on quality scores. This\n" \
    "               has the same effect as the BWA parameter -q:\n" \
    "               BWA trims a read down to argmax_x{sum_{i=x+1}^l(INT-q_i)} \n" \
    "               if q_l<INT where l is the original read length (default 10).\n" \
    "--oldillumina  Specifies that quality score are given in old Phred\n" \
    "               format (Illumina 1.3+, Illumina 1.5+) where quality\n" \
    "               scores are given as Phred + 64 rather than Phred + 33\n" \
    "--rg-sample <STRING>  Use this in the read group SM tag\n" \
    "--rg-lib <STRING>     Use this in the read group LB tag\n" \
    "--multi        Report reads mapping to multiple genomic locations.\n" \
    "               Alternate alignments given in XA tag\n" \
    "--noweb        Do not report any user information and paramters to Amazon S3.\n" \
    "\n\nAdvanced options - general:\n" \
    "-p,--threads <INT>         number of threads (default: 1)\n" \
    "--min-read-length <INT>    minimum number of nucleotides for a\n" \
    "                           read to be processed.\n" \
    "                           (default: 45)\n" \
    "--max-read-length <INT>    maximum number of nucleotides for a\n" \
    "                           read to be processed. (default: 1024)\n" \
    "\n\nAdvanced options - detection:\n"                               \
    "--fft-window-size <INT>    size of fft window (default: 24)\n" \
    "--fft-window-step <INT>    step size of sliding window\n" \
    "                           (default: 12)\n" \
    "--entropy-threshold <FLOAT> threshold score to call a window periodic\n"\
    "                            (defualt: 0.45)\n" \
    "--minperiod <INT>          minimum period to attempt to detect\n"\
    "                           must be >= 2 (default: 2)\n" \
    "--maxperiod <INT>          maximum period to attempt to detect\n"\
    "                           must be <= 6 (default: 6)\n" \
    "--minflank <INT>           minimum length of flanking region to\n" \
    "                           try to align (default: 10)\n" \
    "--maxflank <INT>           length to trim the ends of flanking\n" \
    "                           regions to if they exceed that length\n" \
    "                           (default: 100)\n" \
    "\n\nAdvanced options - alignment:\n" \
    "--max-diff-ref <INT>       maximum difference in length from\n" \
    "                           the reference sequence to report\n" \
    "                           (default: 50bp)\n" \
    "--extend <INT>             Number of bp the reference was extended\n" \
    "                           when building the index.\n" \
    "                           Must be same as --extend parameter used \n" \
    "                           to run lobstr_index.py\n" \
    "                           alignment (default: 1000)\n" \
    "--mapq <INT>               maximum allowed mapq score (default: 100)\n" \
    "                           calculated as the sum of qualities at base \n" \
    "                           mismatches.\n" \
    "-u                         require length difference to be a\n" \
    "                           multiple of the repeat unit\n" \
    "-m,--mismatch <int>        edit distance allowed during alignment\n" \
    "                           of each flanking region (default: -1)\n" \
    "-g <int>                   maximum number of gap opens allowed\n" \
    "                           in each flanking region (default: 1)\n" \
    "-e <int>                   maximum number of gap extensions\n" \
    "                           allowed in each flanking region\n" \
    "                           (default: 1)\n" \
    "-r <float>                 edit distance allowed during alignment\n" \
    "                           of each flanking region (ignored if -m\n" \
    "                           is set) (default: 0.01)\n" \
    "\n\nAdvanced options - Amazon Web Services:\n" \
    "--use-s3 <bucket>          Files are read from this s3 bucket\n" \
    "                           WARNING s3 mode DELETES FILES after processing\n" \
    "                           DO NOT USE this option unless you are pulling\n" \
    "                           files from Amazon S3!\n" \
    "--s3config <file>          s3cmd configuration file (created by\n" \
    "                           s3cmd --configure)\n" \
    "This program takes in raw reads, detects and aligns reads\n" \
    "containing microsatellites, and genotypes STR locations.\n\n";
  cerr << help;
  exit(1);
}

/*
 * parse the command line options
 */
void parse_commandline_options(int argc, char* argv[]) {
  enum LONG_OPTIONS {
    OPT_FILES,
    OPT_PAIR1,
    OPT_PAIR2,
    OPT_GZIP,
    OPT_TABLE,
    OPT_GENOME,
    OPT_OUTPUT,
    OPT_HELP,
    OPT_VERBOSE,
    OPT_DEBUG,
    OPT_ALIGN_DEBUG,
    OPT_FASTQ,
    OPT_BAM,
    OPT_BAMPAIR,
    OPT_THREADS,
    OPT_MISMATCH,
    OPT_NOWEB,
    OPT_RMDUP,
    OPT_FFT_WINDOW_SIZE,
    OPT_FFT_WINDOW_STEP,
    OPT_LOBE_THRESHOLD,
    OPT_EXTEND,
    OPT_MAX_PERIOD,
    OPT_MIN_PERIOD,
    OPT_MIN_FLANK_LEN,
    OPT_MAX_FLANK_LEN,
    OPT_MAX_DIFF_REF,
    OPT_MULTI,
    OPT_FFTW_DEBUG,
    OPT_LOBE_DEBUG,
    OPT_MIN_READ_LENGTH,
    OPT_MAX_READ_LENGTH,
    OPT_ERROR_RATE,
    OPT_MIN_COVERAGE,
    OPT_WHY_NOT,
    OPT_ENTROPY_THRESHOLD,
    OPT_DEBUG_ENTROPY,
    OPT_PROFILE,
    OPT_INDEX,
    OPT_GAP_OPEN,
    OPT_GAP_EXTEND,
    OPT_FPR,
    OPT_EXTEND_FLANK,
    OPT_SW,
    OPT_DEBUGADJUST,
    OPT_ORIG_READ,
    OPT_MAPQ,
    OPT_BWAQ,
    OPT_OLDILLUMINA,
    OPT_CHECKNEXTBEST,
    OPT_USES3,
    OPT_S3CONFIG,
    OPT_S3DEBUG,
    OPT_RG_SAMPLE,
    OPT_RG_LIB,
    OPT_VERSION,
  };

  int ch;
  int option_index = 0;

  static struct option long_options[] = {
    {"files", 1, 0, OPT_FILES},
    {"p1", 1, 0, OPT_PAIR1},
    {"p2", 1, 0, OPT_PAIR2},
    {"gzip", 0, 0, OPT_GZIP},
    {"table", 1, 0, OPT_TABLE},
    {"genome", 1, 0, OPT_GENOME},
    {"out", 1, 0, OPT_OUTPUT},
    {"threads", 1, 0, OPT_THREADS},
    {"noweb", 1, 0, OPT_NOWEB},
    {"mismatch", 1, 0, OPT_MISMATCH},
    {"fft-window-size", 1, 0, OPT_FFT_WINDOW_SIZE},
    {"fft-window-step", 1, 0, OPT_FFT_WINDOW_STEP},
    {"lobe-threshold", 1, 0, OPT_LOBE_THRESHOLD},
    {"extend", 1, 0, OPT_EXTEND},
    {"maxperiod", 1, 0, OPT_MAX_PERIOD},
    {"minperiod", 1, 0, OPT_MIN_PERIOD},
    {"minflank", 1, 0, OPT_MIN_FLANK_LEN},
    {"maxflank", 1, 0, OPT_MAX_FLANK_LEN},
    {"max-diff-ref", 1, 0, OPT_MAX_DIFF_REF},
    {"multi", 0, 0, OPT_MULTI},
    {"help", 0, 0, OPT_HELP},
    {"verbose", 0, 0, OPT_VERBOSE},
    {"debug", 0, 0, OPT_DEBUG},
    {"fastq", 0, 0, OPT_FASTQ},
    {"bam", 0, 0, OPT_BAM},
    {"bampair", 0, 0, OPT_BAMPAIR},
    {"fftw-debug", 0, 0, OPT_FFTW_DEBUG},
    {"lobe-debug", 0, 0, OPT_LOBE_DEBUG},
    {"align-debug", 0, 0, OPT_ALIGN_DEBUG},
    {"min-read-length", 1, 0, OPT_MIN_READ_LENGTH},
    {"max-read-length", 1, 0, OPT_MAX_READ_LENGTH},
    {"min-coverage", 1, 0, OPT_MIN_COVERAGE},
    {"entropy-threshold", 1, 0, OPT_ENTROPY_THRESHOLD},
    {"entropy-debug", 0, 0, OPT_DEBUG_ENTROPY},
    {"profile", 0, 0, OPT_PROFILE},
    {"index-prefix", 1, 0, OPT_INDEX},
    {"extend-flank", 1, 0, OPT_EXTEND_FLANK},
    {"nw-score", 1, 0, OPT_SW},
    {"mapq", 1, 0, OPT_MAPQ},
    {"bwaq", 1, 0, OPT_BWAQ},
    {"oldillumina", 0, 0, OPT_OLDILLUMINA},
    {"debug-adjust", 0, 0, OPT_DEBUGADJUST},
    {"why-not", 0, 0, OPT_WHY_NOT},
    {"orig", 0, 0, OPT_ORIG_READ},
    {"nextbest", 0, 0, OPT_CHECKNEXTBEST},
    {"use-s3", 1, 0, OPT_USES3},
    {"s3config", 1, 0, OPT_S3CONFIG},
    {"s3debug", 0, 0, OPT_S3DEBUG},
    {"rg-sample", 1, 0, OPT_RG_SAMPLE},
    {"rg-lib", 1, 0, OPT_RG_LIB},
    {"version", 0, 0, OPT_VERSION},
    {NULL, no_argument, NULL, 0},
  };
  program = LOBSTR;
  ch = getopt_long(argc, argv, "hvqp:f:t:g:o:m:s:d:e:g:r:u?",
                   long_options, &option_index);
  while (ch != -1) {
    switch (ch) {
    case 'u':
      unit++;
      AddOption("unit", "", false, &user_defined_arguments);
      break;
    case OPT_PROFILE:
      profile++;
      break;
    case OPT_WHY_NOT:
      why_not_debug++;
      break;
    case OPT_ALIGN_DEBUG:
      align_debug++;
      break;
    case OPT_DEBUG:
      debug++;
      break;
    case OPT_DEBUGADJUST:
      debug_adjust++;
      break;
    case 'v':
    case OPT_VERBOSE:
      my_verbose++;
    break;
    case 'h':
    case OPT_HELP:
      show_help();
      exit(1);
    case 'q':
    case OPT_FASTQ:
      input_type = INPUT_FASTQ;
      fastq++;
      AddOption("fastq", "", false, &user_defined_arguments);
      break;
    case OPT_BAM:
      input_type = INPUT_BAM;
      bam++;
      AddOption("bam", "", false, &user_defined_arguments);
      break;
    case OPT_BAMPAIR:
      paired = true;
      input_type = INPUT_BAM;
      bam++;
      AddOption("bampair", "", false, &user_defined_arguments);
      break;
    case 'p':
    case OPT_THREADS:
      threads = atoi(optarg);
      AddOption("threads", string(optarg), true, &user_defined_arguments);
      if (threads <= 0) {
        PrintMessageDieOnError("Invalid number of threads", ERROR);
      }
      if (threads > 1) {
        PrintMessageDieOnError("Multithreading on very small files may fail to produce BAM output", WARNING);
      }
      break;
    case OPT_NOWEB:
      noweb++;
      AddOption("noweb", "", false, &user_defined_arguments);
      break;
    case 'f':
    case OPT_FILES:
      input_files_string = optarg;
      AddOption("files", string(optarg), true, &user_defined_arguments);
      break;
    case OPT_PAIR1:
      input_files_string_p1 = optarg;
      paired = true;
      AddOption("files1", string(optarg), true, &user_defined_arguments);
      break;
    case OPT_PAIR2:
      input_files_string_p2 = optarg;
      paired = true;
      AddOption("files2", string(optarg), true, &user_defined_arguments);
      break;
    case OPT_GZIP:
      user_defined_arguments += "input_gzipped;";
      gzip++;
      AddOption("gzip", "", false, &user_defined_arguments);
      break;
    case 'o':
    case OPT_OUTPUT:
      output_prefix = string(optarg);
      AddOption("out", string(optarg), true, &user_defined_arguments);
      break;
    case 'm':
    case OPT_MISMATCH:
      allowed_mismatches = atoi(optarg);
      if (allowed_mismatches < 0) {
        PrintMessageDieOnError("Invalid number of mismatches", ERROR);
      }
      AddOption("m", string(optarg), true, &user_defined_arguments);
      break;
    case 'b':
    case OPT_FFT_WINDOW_SIZE:
      fft_window_size = atoi(optarg);
      AddOption("fft-window-size", string(optarg), true, &user_defined_arguments);
      break;
    case OPT_FFT_WINDOW_STEP:
      fft_window_step = atoi(optarg);
      AddOption("fft-window-step", string(optarg), true, &user_defined_arguments);
      break;
    case OPT_EXTEND:
      extend = atoi(optarg);
      if (extend <= 0) {
        PrintMessageDieOnError("Invalid extension length", ERROR);
      }
      AddOption("extend", string(optarg), true, &user_defined_arguments);
      break;
    case OPT_MIN_PERIOD:
      min_period = atoi(optarg);
      if (min_period <= 0) {
        PrintMessageDieOnError("Invalid min period", ERROR);
      }
      AddOption("minperiod", string(optarg), true, &user_defined_arguments);
      break;
    case OPT_MAX_PERIOD:
      max_period = atoi(optarg);
      if (max_period <= 0) {
        PrintMessageDieOnError("Invalid max period", ERROR);
      }
      AddOption("maxperiod", string(optarg), true, &user_defined_arguments);
      break;
    case OPT_MAX_FLANK_LEN:
      max_flank_len = atoi(optarg);
      if (max_flank_len <= 0) {
        PrintMessageDieOnError("Invalid max flank length", ERROR);
      }
      AddOption("maxflank", string(optarg), true, &user_defined_arguments);
      break;
    case OPT_MIN_FLANK_LEN:
      min_flank_len = atoi(optarg);
      if (min_flank_len <= 0) {
        PrintMessageDieOnError("Invalid min flank length", ERROR);
      }
      AddOption("minflank", string(optarg), true, &user_defined_arguments);
      break;
    case OPT_MAX_DIFF_REF:
      max_diff_ref = atoi(optarg);
      if (max_diff_ref <=0 ) {
        PrintMessageDieOnError("Invalid max diff ref", ERROR);
      }
      AddOption("max-diff-ref", string(optarg), true, &user_defined_arguments);
      break;
    case OPT_MULTI:
      allow_multi_mappers++;
      AddOption("multi", "", false, &user_defined_arguments);
      break;
    case OPT_FFTW_DEBUG:
      fftw_debug = true;
      break;
    case OPT_LOBE_DEBUG:
      lobe_debug = true;
      break;
    case OPT_MIN_READ_LENGTH:
      min_read_length = atoi(optarg);
      AddOption("min-read-length", string(optarg), true, &user_defined_arguments);
      break;
    case OPT_MAX_READ_LENGTH:
      max_read_length = atoi(optarg);
      AddOption("max-read-length", string(optarg), true, &user_defined_arguments);
      break;
    case OPT_ENTROPY_THRESHOLD:
      entropy_threshold = atof(optarg);
      AddOption("entropy-threshold", string(optarg), true, &user_defined_arguments);
      break;
    case OPT_DEBUG_ENTROPY:
      entropy_debug++;
      break;
    case OPT_INDEX:
      index_prefix = string(optarg);
      AddOption("index-prefix", string(optarg), true, &user_defined_arguments);
      break;
    case OPT_GAP_OPEN:
    case 'g':
      gap_open = atoi(optarg);
      AddOption("gap-open", string(optarg), true, &user_defined_arguments);
      break;
    case OPT_GAP_EXTEND:
    case 'e':
      gap_extend = atoi(optarg);
      AddOption("gap-extend", string(optarg), true, &user_defined_arguments);
      break;
    case OPT_FPR:
    case 'r':
      fpr = atof(optarg);
      AddOption("fpr", string(optarg), true, &user_defined_arguments);
      break;
    case OPT_EXTEND_FLANK:
      extend_flank = atoi(optarg);
      AddOption("extend-flank", string(optarg), true, &user_defined_arguments);
      break;
    case OPT_SW:
      min_sw_score = atoi(optarg);
      AddOption("min-sw-score", string(optarg), true, &user_defined_arguments);
      break;
    case OPT_MAPQ:
      max_mapq = atoi(optarg);
      AddOption("mapq", string(optarg), true, &user_defined_arguments);
      break;
    case OPT_BWAQ:
      QUAL_CUTOFF = atoi(optarg);
      AddOption("bwaq", string(optarg), true, &user_defined_arguments);
      break;
    case OPT_OLDILLUMINA:
      QUALITY_CONSTANT = 64;
      AddOption("oldillumina", "", false, &user_defined_arguments);
      break;
    case OPT_ORIG_READ:
      include_orig_read_start++;
      user_defined_arguments += "orig-read-start=True;";
      break;
    case OPT_CHECKNEXTBEST:
      check_next_best++;
      AddOption("check-next-best", "", false, &user_defined_arguments);
      break;
    case OPT_USES3:
      using_s3++;
      s3bucket = string(optarg);
      AddOption("use-s3", string(optarg), true, &user_defined_arguments);
      break;
    case OPT_S3CONFIG:
      s3cmd_configfile = string(optarg);
      using_s3++;
      AddOption("s3cmdconfig", string(optarg), true, &user_defined_arguments);
      break;
    case OPT_S3DEBUG:
      s3debug++;
      user_defined_arguments += "s3debug;";
      break;
    case OPT_RG_SAMPLE:
      read_group_sample = string(optarg);
      AddOption("rg-sample", string(optarg), true, &user_defined_arguments);
      break;
    case OPT_RG_LIB:
      read_group_library = string(optarg);
      AddOption("rg-lib", string(optarg), true, &user_defined_arguments);
      break;
    case OPT_VERSION:
      cerr << _GIT_VERSION << endl;
      exit(0);
    case '?':
      show_help();
      exit(1);
      break;
    default:
      show_help();
      exit(1);
    }
    ch = getopt_long(argc, argv, "hvqp:f:t:g:o:m:s:d:e:g:r:u?",
                     long_options, &option_index);
  }

  // any arguments left over are extra
  if (optind < argc) {
    PrintMessageDieOnError("Unnecessary leftover arguments", ERROR);
  }
  // make sure arguments make sense
  if (fft_window_step > fft_window_size) {
    PrintMessageDieOnError("fft_window_step must be <=fft_window_size", ERROR);
  }
  if (min_period > max_period) {
    PrintMessageDieOnError("min_period must be <=max_period", ERROR);
  }
  if (min_flank_len > max_flank_len) {
    PrintMessageDieOnError("min_flank_len must be <= max_flank_len", ERROR);
  }
  if (min_period < 1 || max_period > 6) {
    PrintMessageDieOnError("lobSTR can currently only profile STRs for periods 1 thorugh 6", ERROR);
  }
  // check that we have the mandatory parameters
  if ((((!paired || bam) && input_files_string.empty()) ||
       (paired && !bam && (input_files_string_p1.empty() ||
                           input_files_string_p2.empty())))||
      output_prefix.empty() || index_prefix.empty()) {
    PrintMessageDieOnError("Required arguments are missing", ERROR);
  }
  if (gzip && bam) {
    PrintMessageDieOnError("Gzip option not compatible with bam input", ERROR);
  }
  if (using_s3 && s3cmd_configfile.empty()) {
    PrintMessageDieOnError("Must supply an s3cmd config file", ERROR);
  }
}

void LoadReference(const std::string& repseq) {
  // Load BWT index
  string prefix = index_prefix;
  prefix += repseq;
  prefix += ".fa";

  string bwt_str = prefix+".bwt";
  bwt_t *bwt_forward, *bwt_reverse;
  bwt_forward = bwt_restore_bwt(bwt_str.c_str());

  string rbwt_str = prefix+".rbwt";
  bwt_reverse = bwt_restore_bwt(rbwt_str.c_str());

  string sa_str = prefix+".sa";
  bwt_restore_sa(sa_str.c_str(), bwt_forward);

  string rsa_str = prefix+".rsa";
  bwt_restore_sa(rsa_str.c_str(), bwt_reverse);

  BWT bwt_ref;
  bwt_ref.bwt[0] = bwt_forward;
  bwt_ref.bwt[1] = bwt_reverse;
  bwt_references.insert(pair<string, BWT>(repseq, bwt_ref));

  // Load BNT annotations
  BNT bnt;
  bntseq_t *bns;
  bns = bns_restore(prefix.c_str());
  bnt.bns = bns;
  bnt_annotations.insert(pair<string, BNT>(repseq, bnt));
}

void DestroyReferences() {
  for (map<string, BWT>::iterator it = bwt_references.begin();
       it != bwt_references.end(); ++it) {
    bwt_destroy(it->second.bwt[0]);
    bwt_destroy(it->second.bwt[1]);
  }
  for (map<string, BNT>::iterator it = bnt_annotations.begin();
       it != bnt_annotations.end(); ++it) {
    bns_destroy(it->second.bns);
  }
}



/*
 * process read in single thread
 */
void single_thread_process_loop(const vector<string>& files1,
                                const vector<string>& files2) {
  ReadPair read_pair;
  SamFileWriter samWriter(output_prefix + ".aligned.bam", chrom_sizes);
  STRDetector *pDetector = new STRDetector();
  BWAReadAligner *pAligner = new BWAReadAligner(&bwt_references,
                                                &bnt_annotations,
                                                &ref_sequences, opts);
  std::string file1;
  std::string file2;
  for (size_t i = 0; i < files1.size(); i++) {
    file1 = files1.at(i);
    if (paired && !bam) {
      file2 = files2.at(i);
      PrintMessageDieOnError("Processing files " + file1 + " and " + file2, PROGRESS);
      if (using_s3) {
        const std::string s3cmd1 = GenerateS3Command(s3bucket,
                                                     file1,
                                                     s3cmd_configfile);
        const std::string s3cmd2 = GenerateS3Command(s3bucket,
                                                     file2,
                                                     s3cmd_configfile);
        if (s3debug) {
          PrintMessageDieOnError("S3 debug: " + s3cmd1, PROGRESS);
          PrintMessageDieOnError("S3 debug: " + s3cmd2, PROGRESS);
        } else {
          if (system(s3cmd1.c_str()) != 0) {
            PrintMessageDieOnError("Problem fetching file1 from S3", ERROR);
          }
          if (system(s3cmd2.c_str()) != 0) {
            PrintMessageDieOnError("Problem fetching file2 from S3", ERROR);
          }
        }
        file1 = "/mnt/lobstr/"+file1;
        file2 = "/mnt/lobstr/"+file2;
      }
      if (!(fexists(file1.c_str()) && fexists(file2.c_str()))) {
        PrintMessageDieOnError("File " + file1 + " or " + file2 + " does not exist", WARNING);
        continue;
      }
    } else {
      PrintMessageDieOnError("Processing file " + file1, PROGRESS);
      if (using_s3) {
        const std::string s3cmd = GenerateS3Command(s3bucket,
                                                    file1,
                                                    s3cmd_configfile);
        if (s3debug) {
          PrintMessageDieOnError("S3 debug: " + s3cmd, PROGRESS);
        } else {
          if (system(s3cmd.c_str()) != 0) {
            PrintMessageDieOnError("Problem fetching file from S3", ERROR);
          }
        }
        file1 = "/mnt/lobstr/"+file1;
      }
      if (!fexists(file1.c_str())) {
        PrintMessageDieOnError("File " + file1 + " does not exist", WARNING);
        continue;
      }
    }
    IFileReader* pReader = create_file_reader(file1, file2);
    int aligned = false;
    int num_reads_processed = 0;
    std::string repseq = "";
    while (pReader->GetNextRecord(&read_pair)) {
      aligned = false;
      num_reads_processed += 1;
      if (num_reads_processed % READPROGRESS == 0) {
        stringstream msg;
        msg << "Processed " << num_reads_processed << " reads";
        PrintMessageDieOnError(msg.str(), PROGRESS);
      }
      read_pair.read_count = num_reads_processed;
      // reset fields
      read_pair.reads.at(0).repseq = "";
      read_pair.reads.at(0).ms_repeat_best_period = 0;
      read_pair.reads.at(0).ms_repeat_next_best_period = 0;
      if (read_pair.reads.at(0).paired) {
        read_pair.reads.at(1).repseq = "";
        read_pair.reads.at(1).ms_repeat_best_period = 0;
        read_pair.reads.at(1).ms_repeat_next_best_period = 0;
      }

      // Check read length
      if (!(read_pair.reads.at(0).nucleotides.length() >= min_read_length) &&
          (read_pair.reads.at(0).nucleotides.length() <= max_read_length)) {
        continue;
      }
      if (read_pair.reads.at(0).paired) {
        if (!(read_pair.reads.at(1).nucleotides.length() >= min_read_length) &&
            (read_pair.reads.at(1).nucleotides.length() <= max_read_length)) {
          continue;
        }
      }
      bases += read_pair.reads.at(0).nucleotides.length();
      if (read_pair.reads.at(0).paired) bases += read_pair.reads.at(1).nucleotides.length();

      // STEP 1: Sensing
      string det_err, det_messages;
      if (!pDetector->ProcessReadPair(&read_pair, &det_err, &det_messages)) {
        if (debug) {
          PrintMessageDieOnError(GetReadDebug(read_pair, det_err, det_messages, "NA", "NA") + " (detection-fail)", DEBUG);
        }
        continue;
      }
      // STEP 2: Alignment
      string aln_err, aln_messages;
      if (pAligner->ProcessReadPair(&read_pair, &aln_err, &aln_messages)) {
        aligned = true;
        if (debug) { // if aligned, what was the repseq we aligned to
          PrintMessageDieOnError(GetReadDebug(read_pair, det_err, det_messages, aln_err, aln_messages)+ " (aligned-round-1)", DEBUG);
        }
      } else {
        read_pair.read1_passed_detection = false;
        read_pair.read2_passed_detection = false;
        // Try second best period for each read
        if (read_pair.reads.at(0).ms_repeat_next_best_period != 0) {
          read_pair.reads.at(0).ms_repeat_best_period =
            read_pair.reads.at(0).ms_repeat_next_best_period;
          string err, second_best_repseq;
          if (getMSSeq(read_pair.reads.at(0).detected_ms_region_nuc,
                       read_pair.reads.at(0).ms_repeat_best_period, &repseq, &second_best_repseq, &err)) {
            read_pair.reads.at(0).repseq = repseq;
            read_pair.read1_passed_detection = true;
          }
        }
        if (read_pair.reads.at(0).paired) {
          if (read_pair.reads.at(1).ms_repeat_next_best_period != 0) {
            read_pair.reads.at(1).ms_repeat_best_period =
              read_pair.reads.at(1).ms_repeat_next_best_period;
            string err, second_best_repseq;
            if (getMSSeq(read_pair.reads.at(1).detected_ms_region_nuc,
                         read_pair.reads.at(1).ms_repeat_best_period,
                         &repseq, &second_best_repseq, &err)) {
              read_pair.reads.at(1).repseq = repseq;
              read_pair.read2_passed_detection = true;
            }
          }
        }
        if (read_pair.read1_passed_detection ||
            read_pair.read2_passed_detection) {
          if (pAligner->ProcessReadPair(&read_pair, &aln_err, &aln_messages)) {
            aligned = true;
            if (debug) { // if aligned, what was the repseq we aligned to
              PrintMessageDieOnError(GetReadDebug(read_pair, det_err, det_messages, aln_err, aln_messages)+ " (aligned-round-2)", DEBUG);
            }
          }
        }
      }
      if (aligned) {
        samWriter.WriteRecord(read_pair);
      } else {
        if (debug) { // if didn't align, print this
          PrintMessageDieOnError(GetReadDebug(read_pair, det_err, det_messages, aln_err, aln_messages)+ " (not-aligned)", DEBUG);
        }
      }
    }
    delete pReader;
    stringstream msg;
    msg << "Processed " << num_reads_processed << " reads";
    PrintMessageDieOnError(msg.str(), PROGRESS);
    if (using_s3) {
      string rmcmd = "rm " + file1;
      if (paired && !bam) {
        rmcmd += "; rm ";
        rmcmd += file2;
      }
      if (s3debug) {
        PrintMessageDieOnError("S3 debug: " + rmcmd, PROGRESS);
      } else {
        if (system(rmcmd.c_str()) != 0) {
          PrintMessageDieOnError("Problem deleting file", ERROR);
        }
      }
    }
  }
  delete pDetector;
  delete pAligner;
}


void* satellite_process_consumer_thread(void *arg) {
  MultithreadData *pMT_DATA = reinterpret_cast<MultithreadData*>(arg);
  STRDetector *pDetector = new STRDetector();
  BWAReadAligner *pAligner = new BWAReadAligner(&bwt_references,
                                                &bnt_annotations,
                                                &ref_sequences, opts);
  int aligned = false;
  while (1) {
    aligned = false;
    std::string repseq;
    ReadPair* pReadRecord = pMT_DATA->get_new_input();
    if (!(pReadRecord->reads.at(0).nucleotides.length() >= min_read_length)
        && (pReadRecord->reads.at(0).nucleotides.length() <= max_read_length)) {
      delete pReadRecord;
      pMT_DATA->increment_output_counter();
      continue;
    }
    if (pReadRecord->reads.at(0).paired) {
      if (!(pReadRecord->reads.at(1).nucleotides.length() >= min_read_length) &&
          (pReadRecord->reads.at(1).nucleotides.length() <= max_read_length)) {
        delete pReadRecord;
        pMT_DATA->increment_output_counter();
        continue;
      }
    }
    bases += pReadRecord->reads.at(0).nucleotides.length();
    if (pReadRecord->reads.at(0).paired) bases += pReadRecord->reads.at(1).nucleotides.length();

    // Reset fields
    pReadRecord->reads.at(0).repseq = "";
    pReadRecord->reads.at(0).ms_repeat_best_period = 0;
    pReadRecord->reads.at(0).ms_repeat_next_best_period = 0;
    if (pReadRecord->reads.at(0).paired) {
      pReadRecord->reads.at(1).repseq = "";
      pReadRecord->reads.at(1).ms_repeat_best_period = 0;
      pReadRecord->reads.at(1).ms_repeat_next_best_period = 0;
    }

    // STEP 1: Sensing
    string err, messages;
    if (!pDetector->ProcessReadPair(pReadRecord, &err, &messages)) {
      pMT_DATA->increment_output_counter();
      delete pReadRecord;
      continue;
    }

    // STEP 2: Alignment
    if (pAligner->ProcessReadPair(pReadRecord, &err, &messages)) {
      aligned = true;
    } else {
      pReadRecord->read1_passed_detection = false;
      pReadRecord->read2_passed_detection = false;
      // Try second best period for each read
      if (pReadRecord->reads.at(0).ms_repeat_next_best_period != 0) {
        pReadRecord->reads.at(0).ms_repeat_best_period =
          pReadRecord->reads.at(0).ms_repeat_next_best_period;
        string err, second_best_repseq;
        if (getMSSeq(pReadRecord->reads.at(0).detected_ms_region_nuc,
                     pReadRecord->reads.at(0).ms_repeat_best_period, &repseq, &second_best_repseq, &err)) {
          pReadRecord->reads.at(0).repseq = repseq;
          pReadRecord->read1_passed_detection = true;
        }
      }
      if (pReadRecord->reads.at(0).paired) {
        if (pReadRecord->reads.at(1).ms_repeat_next_best_period != 0) {
          pReadRecord->reads.at(1).ms_repeat_best_period =
            pReadRecord->reads.at(1).ms_repeat_next_best_period;
          string err, second_best_repseq;
          if (getMSSeq(pReadRecord->reads.at(1).detected_ms_region_nuc,
                       pReadRecord->reads.at(1).ms_repeat_best_period,
                       &repseq, &second_best_repseq, &err)) {
            pReadRecord->reads.at(1).repseq = repseq;
            pReadRecord->read2_passed_detection = true;
          }
        }
      }
      if (pReadRecord->read1_passed_detection ||
          pReadRecord->read2_passed_detection) {
        if (pAligner->ProcessReadPair(pReadRecord, &err, &messages)) {
          aligned = true;
        }
      }
    }
    if (aligned) {
      pMT_DATA->post_new_output_read(pReadRecord);
    } else {
      // Don't pass this read on to the output thread,
      // discard it here
      delete pReadRecord;
      pMT_DATA->increment_output_counter();
    }
  }
  pthread_exit(reinterpret_cast<void*>(arg));
}

void* output_writer_thread(void *arg) {
  MultithreadData *pMT_DATA = reinterpret_cast<MultithreadData*>(arg);
  SamFileWriter samWriter(output_prefix + ".aligned.bam", chrom_sizes);
  while (1) {
    ReadPair *pReadRecord = pMT_DATA->get_new_output();
    samWriter.WriteRecord(*pReadRecord);
    delete pReadRecord;
    pMT_DATA->increment_output_counter();
  }
  delete pMT_DATA;
  pthread_exit(reinterpret_cast<void*>(arg));
}

void multi_thread_process_loop(vector<string> files1,
                               vector<string> files2) {
  MultithreadData mtdata(threads);
  list<pthread_t> satellite_threads;
  pthread_t writer_thread;
  if (files1.size() == 0) return;
  for (size_t i = 0; i < threads; ++i) {
    pthread_t id;
    if (pthread_create(&id, NULL, satellite_process_consumer_thread,
                       reinterpret_cast<void*>(&mtdata))) {
      PrintMessageDieOnError("Failed to create threads", ERROR);
    }
    satellite_threads.push_back(id);
  }

  if (pthread_create(&writer_thread, NULL, output_writer_thread,
                     reinterpret_cast<void*>(&mtdata))) {
    PrintMessageDieOnError("Failed to create output writer threads", ERROR);
  }

  size_t counter = 1;
  std::string file1;
  std::string file2;
  for (size_t i = 0; i < files1.size(); i++) {
    file1 = files1.at(i);
    if (paired && !bam) {
      file2 = files2.at(i);
      PrintMessageDieOnError("Processing files " + file1 + " and " + file2, PROGRESS);
      if (using_s3) {
        const std::string s3cmd1 = GenerateS3Command(s3bucket,
                                                     file1,
                                                     s3cmd_configfile);
        const std::string s3cmd2 = GenerateS3Command(s3bucket,
                                                     file2,
                                                     s3cmd_configfile);
        if (s3debug) {
          PrintMessageDieOnError("S3 debug: " + s3cmd1, PROGRESS);
          PrintMessageDieOnError("S3 debug: " + s3cmd2, PROGRESS);
        } else {
          if (system(s3cmd1.c_str()) != 0) {
            PrintMessageDieOnError("Problem fetching file1 from S3", ERROR);
          }
          if (system(s3cmd2.c_str()) != 0) {
            PrintMessageDieOnError("Problem fetching file2 from S3", ERROR);
          }
        }
        file1 = "/mnt/lobstr/"+file1;
        file2 = "/mnt/lobstr/"+file2;
      }
      if (!(fexists(file1.c_str()) && fexists(file2.c_str()))) {
        PrintMessageDieOnError("File " + file1 + " or " + file2 + " does not exist", WARNING);
        continue;
      }
    } else {
      PrintMessageDieOnError("Processing file " + file1, PROGRESS);
      if (using_s3) {
        const std::string s3cmd = GenerateS3Command(s3bucket,
                                                    file1,
                                                    s3cmd_configfile);
        if (s3debug) {
          PrintMessageDieOnError("S3 debug: " + s3cmd, PROGRESS);
        } else {
          if (system(s3cmd.c_str()) != 0) {
            PrintMessageDieOnError("Problem fetching file from S3", ERROR);
          }
        }
        file1 = "/mnt/lobstr/"+file1;
        file2 = "/mnt/lobstr/"+file2;
      }
      if (!fexists(file1.c_str())) {
        PrintMessageDieOnError("File " + file1 + " or " + file2 + " does not exist", WARNING);
        continue;
      }
    }
    IFileReader *pReader = create_file_reader(file1, file2);
    do {
      ReadPair *pRecord = new ReadPair;
      pRecord->read_count = counter;
      if (counter % READPROGRESS == 0) {
        stringstream msg;
        msg << "Processed " << counter << " reads";
        PrintMessageDieOnError(msg.str(), PROGRESS);
      }
      if (!pReader->GetNextRecord(pRecord))
        break;  // no more reads
      counter++;
      mtdata.increment_input_counter();
      mtdata.post_new_input_read(pRecord);
      pRecord = NULL;  // the consumers will take it from here, and free it
    } while (1);
    delete pReader;
    if (using_s3) {
      string rmcmd = "rm " + file1;
      if (paired && !bam) {
        rmcmd += "; rm ";
        rmcmd += file2;
      }
      if (s3debug) {
        PrintMessageDieOnError("S3 debug: " + rmcmd, PROGRESS);
      } else {
        if (system(rmcmd.c_str()) != 0) {
          PrintMessageDieOnError("Problem deleting file", ERROR);
        }
      }
    }
  }
  while (1) {
    sleep(1);  // OMG, the horror...
    if ( mtdata.input_output_counters_equal())
      break;
    sleep(10);
    break;
  }
}

int main(int argc, char* argv[]) {
  PrintLobSTR();
  time_t starttime, endtime;
  time(&starttime);
  parse_commandline_options(argc, argv);
  PrintMessageDieOnError("Getting run info", PROGRESS);
  run_info.Reset();
  run_info.starttime = GetTime();
  if (_GIT_VERSION != NULL) {
    run_info.gitversion = _GIT_VERSION;
  } else {run_info.gitversion = "Not available";}
  if (_MACHTYPE != NULL) {
    run_info.machtype = _MACHTYPE;
  } else {run_info.machtype = "Not available";}
  run_info.params = user_defined_arguments;

  PrintMessageDieOnError("Initializing...", PROGRESS);
  // open file with all names
  TextFileReader tReader(index_prefix+"strdict.txt");
  string line = "";
  while (tReader.GetNextLine(&line)) {
    // set chrom sizes
    if (line.substr(0, 3) == "REF") {
      vector<string> items(0);
      split(line, '\t', items);
      chrom_sizes[items[1]] = atoi(items[2].c_str());
    } else {
      // make sure repeat is valid
      if (count(line, line.at(0)) != line.length() || line.length() == 1) {
        LoadReference(line);
      }
    }
  }

  // Load fasta reference with all STR sequences
  FastaFileReader faReader(index_prefix+"ref.fa");
  MSReadRecord ref_record;
  while (faReader.GetNextRead(&ref_record)) {
    REFSEQ refseq;
    vector<string> items;
    string refstring = ref_record.ID;
    split(refstring, '$', items);
    if (items.size() == 7) {
      refseq.sequence = ref_record.orig_nucleotides;
      refseq.start = atoi(items.at(2).c_str());
      int refid = atoi(items.at(0).c_str());
      ref_sequences.insert(pair<int, REFSEQ>(refid, refseq));
    }
  }

  // set up options
  opts = gap_init_opt();
  opts->max_diff = allowed_mismatches;
  if (allowed_mismatches == -1) {
    opts->fnr = fpr;
  } else {
    opts->fnr = -1.0;
  }
  // only 1 indel per flank
  opts->max_gapo = gap_open;
  opts->max_gape = gap_extend;
  // take the first INT subsequence as seed
  opts->seed_len = 5;

  // get the input files
  if (paired && !bam) {
    boost::split(input_files1, input_files_string_p1, boost::is_any_of(","));
    boost::split(input_files2, input_files_string_p2, boost::is_any_of(","));
    if (!input_files1.size() == input_files2.size()) {
      PrintMessageDieOnError("Different number of files for each pair", ERROR);
    }
  } else {
    boost::split(input_files, input_files_string, boost::is_any_of(","));
  }
  // Initialize fft
  // Create the singleton HammingWindow (before starting any threads)
  HammingWindowGenerator* hamgen =
    HammingWindowGenerator::GetHammingWindowSingleton();
  TukeyWindowGenerator* tukgen =
    TukeyWindowGenerator::GetTukeyWindowSingleton();

  // Initialize global FFTW plans
  FFT_NUC_VECTOR::initialize_fftw_plans();

  // run detection/alignment
  PrintMessageDieOnError("Running detection/alignment...", PROGRESS);
  if (threads == 1) {
    if (paired && !bam) {
      single_thread_process_loop(input_files1, input_files2);
    } else {
      single_thread_process_loop(input_files, vector<string>(0));
    }
  } else {
    if (paired && !bam) {
      multi_thread_process_loop(input_files1, input_files2);
    } else {
      multi_thread_process_loop(input_files, vector<string>(0));
    }
  }
  delete hamgen;
  delete tukgen;
  run_info.endtime = GetTime();
  OutputRunStatistics();
  time(&endtime);
  stringstream msg;
  int seconds_elapsed = difftime(endtime, starttime);
  msg << "Done! " << seconds_elapsed/60/60/24 << ":"
      << (seconds_elapsed/60/60)%24 << ":"
      << (seconds_elapsed/60)%60 << ":"
      << seconds_elapsed%60 << " elapsed";
  PrintMessageDieOnError(msg.str(), PROGRESS);
  return 0;
}
