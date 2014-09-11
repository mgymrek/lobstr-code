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

//Enable the following to write some thread-related debug messages
//#define DEBUG_THREADS

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
#include <sstream>
#include <vector>
#include <utility>

#include "src/BWAReadAligner.h"
#include "src/bwtaln.h"
#include "src/bwase.h"
#include "src/common.h"
#include "src/FastaFileReader.h"
#include "src/FastqFileReader.h"
#include "src/IFileReader.h"
#include "src/MSReadRecord.h"
#include "src/MultithreadData.h"
#include "src/SamFileWriter.h"
#include "src/STRDetector.h"
#include "src/runtime_parameters.h"

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

// Either "reads" or "pairs", depending on what's being processed.
std::string unit_name;

// alignment references, keep global
BNT bnt_annotation;
BWT bwt_reference;
void LoadReference();
void DestroyReference();
gap_opt_t *opts;

void show_help() {
  const char* help =
    "\nlobSTR [OPTIONS] " \
    "    {-f <file1[,file2,...]> | --p1 <file1_1[,file2_1,...]>\n" \
    "    --p2 <file1_2[,file2_1,...]>} --index-prefix <index prefix>\n" \
    "    -o <output prefix> --rg-sample <STRING> --rg-lib <STRING>\n" \
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
    "--rg-sample <STRING>  Use this in the read group SM tag\n" \
    "--rg-lib <STRING>     Use this in the read group LB tag\n" \
    "\n\nOptions:\n" \
    "-h,--help      display this help screen\n" \
    "-v,--verbose   print out useful progress messages\n" \
    "--quiet    don't print anything to stderr or stdout\n" \
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
    "--max-hits-quit-aln <int>  Stop alignment search after int hits found.\n" \
    "                           Default: 1000. Use -1 for no limit.\n" \
    "--min-flank-allow-mismatch <int>  Mininum length of flanking region to allow\n" \
    "                           mismatches. Default: 30.\n" \
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
    OPT_GENOME,
    OPT_OUTPUT,
    OPT_HELP,
    OPT_VERBOSE,
    OPT_QUIET,
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
    OPT_MIN_FLANK_ALLOW_MISMATCH,
    OPT_MAX_HITS_QUIT_ALN,
    OPT_MIN_FLANK_LEN,
    OPT_MAX_FLANK_LEN,
    OPT_MAX_DIFF_REF,
    OPT_MULTI,
    OPT_MIN_READ_LENGTH,
    OPT_MAX_READ_LENGTH,
    OPT_ERROR_RATE,
    OPT_ENTROPY_THRESHOLD,
    OPT_INDEX,
    OPT_GAP_OPEN,
    OPT_GAP_EXTEND,
    OPT_FPR,
    OPT_SW,
    OPT_MAPQ,
    OPT_BWAQ,
    OPT_OLDILLUMINA,
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
    {"genome", 1, 0, OPT_GENOME},
    {"out", 1, 0, OPT_OUTPUT},
    {"threads", 1, 0, OPT_THREADS},
    {"noweb", 0, 0, OPT_NOWEB},
    {"mismatch", 1, 0, OPT_MISMATCH},
    {"fft-window-size", 1, 0, OPT_FFT_WINDOW_SIZE},
    {"fft-window-step", 1, 0, OPT_FFT_WINDOW_STEP},
    {"lobe-threshold", 1, 0, OPT_LOBE_THRESHOLD},
    {"extend", 1, 0, OPT_EXTEND},
    {"minflank", 1, 0, OPT_MIN_FLANK_LEN},
    {"maxflank", 1, 0, OPT_MAX_FLANK_LEN},
    {"max-diff-ref", 1, 0, OPT_MAX_DIFF_REF},
    {"multi", 0, 0, OPT_MULTI},
    {"help", 0, 0, OPT_HELP},
    {"verbose", 0, 0, OPT_VERBOSE},
    {"quiet", 0, 0, OPT_QUIET},
    {"debug", 0, 0, OPT_DEBUG},
    {"fastq", 0, 0, OPT_FASTQ},
    {"bam", 0, 0, OPT_BAM},
    {"bampair", 0, 0, OPT_BAMPAIR},
    {"align-debug", 0, 0, OPT_ALIGN_DEBUG},
    {"min-read-length", 1, 0, OPT_MIN_READ_LENGTH},
    {"max-read-length", 1, 0, OPT_MAX_READ_LENGTH},
    {"min-flank-allow-mismatch", 1, 0, OPT_MIN_FLANK_ALLOW_MISMATCH},
    {"max-hits-quit-aln", 1, 0, OPT_MAX_HITS_QUIT_ALN},
    {"entropy-threshold", 1, 0, OPT_ENTROPY_THRESHOLD},
    {"index-prefix", 1, 0, OPT_INDEX},
    {"nw-score", 1, 0, OPT_SW},
    {"mapq", 1, 0, OPT_MAPQ},
    {"bwaq", 1, 0, OPT_BWAQ},
    {"oldillumina", 0, 0, OPT_OLDILLUMINA},
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
    case OPT_ALIGN_DEBUG:
      align_debug++;
      break;
    case OPT_DEBUG:
      debug++;
      break;
    case 'v':
    case OPT_VERBOSE:
      my_verbose++;
    break;
    case OPT_QUIET:
      quiet++;
      break;
    case 'h':
    case OPT_HELP:
      show_help();
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
    case OPT_MAX_HITS_QUIT_ALN:
      max_hits_quit_aln = atoi(optarg);
      AddOption("max-hits-quit-aln", string(optarg), true, &user_defined_arguments);
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
    default:
      show_help();
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
  if (min_flank_len > max_flank_len) {
    PrintMessageDieOnError("min_flank_len must be <= max_flank_len", ERROR);
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
  if (read_group_sample.empty() || read_group_library.empty()) {
    PrintMessageDieOnError("Must specify --rg-lib and --rg-sample", ERROR);
  }
}

void LoadChromSizes() {
  TextFileReader tReader(index_prefix+"chromsizes.tab");
  string line;
  while (tReader.GetNextLine(&line)) {
    vector<string> items(0);
    split(line, '\t', items);
    if (items.size() != 2) {
      PrintMessageDieOnError("Chromosome sizes file malformed", ERROR);
    }
    chrom_sizes[items[0]] = atoi(items[1].c_str());
  }
}

void LoadReference() {
  // Load BWT index
  string prefix = index_prefix + "ref.fasta";

  string bwt_str = prefix+".bwt";
  bwt_t *bwt_forward, *bwt_reverse;
  bwt_forward = bwt_restore_bwt(bwt_str.c_str());

  string rbwt_str = prefix+".rbwt";
  bwt_reverse = bwt_restore_bwt(rbwt_str.c_str());

  string sa_str = prefix+".sa";
  bwt_restore_sa(sa_str.c_str(), bwt_forward);

  string rsa_str = prefix+".rsa";
  bwt_restore_sa(rsa_str.c_str(), bwt_reverse);

  bwt_reference.bwt[0] = bwt_forward;
  bwt_reference.bwt[1] = bwt_reverse;

  // Load BNT annotations
  bntseq_t *bns;
  bns = bns_restore(prefix.c_str());
  bnt_annotation.bns = bns;

  // Load fasta reference with all STR sequences
  FastaFileReader faReader(index_prefix+"ref.fasta");
  MSReadRecord ref_record;
  while (faReader.GetNextRead(&ref_record)) {
    REFSEQ refseq;
    vector<string> items;
    string refstring = ref_record.ID;
    split(refstring, '$', items);
    if (items.size() == 4) {
      refseq.sequence = ref_record.orig_nucleotides;
      refseq.start = atoi(items.at(2).c_str());
      refseq.chrom = items.at(1);
      int refid = atoi(items.at(0).c_str());
      ref_sequences.insert(pair<int, REFSEQ>(refid, refseq));
    } else {
      PrintMessageDieOnError("Malformed reference fasta", ERROR);
    }
  }

  // Load map
  TextFileReader tReader(index_prefix+"ref_map.tab");
  string line;
  while (tReader.GetNextLine(&line)) {
    vector<string>items(0);
    split(line, '\t', items);
    if (items.size() != 3) {
      PrintMessageDieOnError("Malformed map file", ERROR);
    }
    // Get refid
    int refid = atoi(items.at(0).c_str());
    // Get motifs
    string motif_string = items.at(1);
    vector<string>motifs(0);
    split(motif_string, ';', motifs);
    ref_sequences[refid].motifs = motifs;
    for (size_t i = 0; i < motifs.size(); i++) {
      vector<ReferenceSTR> vec;
      ref_sequences[refid].ref_strs.insert(pair<string, vector<ReferenceSTR> >(motifs.at(i), vec));
    }
    // Get reference STRs
    string regions_string = items.at(2);
    vector<string>regions(0);
    split(regions_string, ';', regions);
    for (size_t i = 0; i < regions.size(); i++) {
      vector<string>region_items(0);
      split(regions.at(i), '_', region_items);
      string motif = region_items.at(3);
      int start = atoi(region_items.at(0).c_str());
      int stop = atoi(region_items.at(1).c_str());
      ReferenceSTR ref_str;
      ref_str.chrom = ref_sequences[refid].chrom;
      ref_str.start = start;
      ref_str.stop = stop;
      ref_sequences[refid].ref_strs[motif].push_back(ref_str);
    }
  }
}

void DestroyReference() {
  bwt_destroy(bwt_reference.bwt[0]);
  bwt_destroy(bwt_reference.bwt[1]);
  bns_destroy(bnt_annotation.bns);
}

/*
 * process read in single thread
 */
void single_thread_process_loop(const vector<string>& files1,
                                const vector<string>& files2) {
  ReadPair read_pair;
  SamFileWriter samWriter(output_prefix + ".aligned.bam", chrom_sizes);
  STRDetector *pDetector = new STRDetector();
  BWAReadAligner *pAligner = new BWAReadAligner(&bwt_reference,
                                                &bnt_annotation,
                                                &ref_sequences, opts);
  std::string file1;
  std::string file2;
  size_t num_reads_processed = 0;
  for (size_t i = 0; i < files1.size(); i++) {
    file1 = files1.at(i);
    if (paired && !bam) {
      file2 = files2.at(i);
      PrintMessageDieOnError("Processing files " + file1 + " and " + file2, PROGRESS);
      if (!(fexists(file1.c_str()) && fexists(file2.c_str()))) {
        PrintMessageDieOnError("File " + file1 + " or " + file2 + " does not exist", WARNING);
        continue;
      }
    } else {
      PrintMessageDieOnError("Processing file " + file1, PROGRESS);
      if (!fexists(file1.c_str())) {
        PrintMessageDieOnError("File " + file1 + " does not exist", WARNING);
        continue;
      }
    }
    IFileReader* pReader = create_file_reader(file1, file2);
    int aligned = false;
    std::string repseq = "";
    while (pReader->GetNextRecord(&read_pair)) {
      aligned = false;
      num_reads_processed += 1;
      if (num_reads_processed % READPROGRESS == 0) {
        stringstream msg;
        msg << "Processed " << num_reads_processed << ' ' << unit_name;
        PrintMessageDieOnError(msg.str(), PROGRESS);
      }
      read_pair.read_count = num_reads_processed;
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
    msg << "Processed " << num_reads_processed << ' ' << unit_name;
    PrintMessageDieOnError(msg.str(), PROGRESS);
  }
  delete pDetector;
  delete pAligner;
  run_info.num_processed_units = num_reads_processed;
}


void* satellite_process_consumer_thread(void *arg) {
  MultithreadData *pMT_DATA = reinterpret_cast<MultithreadData*>(arg);
  STRDetector *pDetector = new STRDetector();
  BWAReadAligner *pAligner = new BWAReadAligner(&bwt_reference,
                                                &bnt_annotation,
                                                &ref_sequences, opts);
  int aligned = false;
#ifdef DEBUG_THREADS
  std::stringstream msg;
  msg << "Alignment thread " << pthread_self() << " started" ;
  PrintMessageDieOnError(msg.str(), PROGRESS);
#endif
  while (1) {
    aligned = false;
    std::string repseq;
    ReadPair* pReadRecord = pMT_DATA->get_new_input();
    if (pReadRecord == NULL) {
#ifdef DEBUG_THREADS
      std::stringstream msg;
      msg << "Alignment thread " << pthread_self() << " completed" ;
      PrintMessageDieOnError(msg.str(), PROGRESS);
#endif
      break;
    }
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
  return NULL;
}

void* output_writer_thread(void *arg) {
  MultithreadData *pMT_DATA = reinterpret_cast<MultithreadData*>(arg);
  SamFileWriter samWriter(output_prefix + ".aligned.bam", chrom_sizes);
#ifdef DEBUG_THREADS
  std::stringstream msg;
  msg << "Writer thread " << pthread_self() << " started (output file ='"
      << output_prefix << ".aligned.bam" << ")" ;
  PrintMessageDieOnError(msg.str(),PROGRESS);
#endif
  while (1) {
    ReadPair *pReadRecord = pMT_DATA->get_new_output();
    if (pReadRecord == NULL) {
#ifdef DEBUG_THREADS
      std::stringstream msg;
      msg << "Writer thread " << pthread_self() << " completed";
      PrintMessageDieOnError(msg.str(),PROGRESS);
#endif
      break;
    }
    samWriter.WriteRecord(*pReadRecord);
    delete pReadRecord;
    pMT_DATA->increment_output_counter();
  }
  return NULL;
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
      if (!(fexists(file1.c_str()) && fexists(file2.c_str()))) {
        PrintMessageDieOnError("File " + file1 + " or " + file2 + " does not exist", WARNING);
        continue;
      }
    } else {
      PrintMessageDieOnError("Processing file " + file1, PROGRESS);
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
        msg << "Processed " << counter << ' ' << unit_name;
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
  }
  run_info.num_processed_units = counter;

#ifdef DEBUG_THREADS
  PrintMessageDieOnError("No more input, waiting for alignment threads completion", PROGRESS);
#endif
  //Send a 'poison pill' to the alignment threads
  for (size_t i = 0; i < threads; ++i)
    mtdata.post_new_input_read(NULL);

  for (list<pthread_t>::const_iterator it = satellite_threads.begin();
          it != satellite_threads.end(); ++it) {
    int i = pthread_join(*it,NULL);
    if (i != 0) {
       stringstream msg;
       msg << "Failed to join alignment thread " << (*it) <<
              "error code = " << i ;
       PrintMessageDieOnError(msg.str(), WARNING);
    }
  }
  //Send a 'poison pill' to the writer thread
#ifdef DEBUG_THREADS
  PrintMessageDieOnError("waiting for writer thread completion", PROGRESS);
#endif
  mtdata.post_new_output_read(NULL);
  int i = pthread_join(writer_thread,NULL);
  if (i != 0) {
    stringstream msg;
    msg << "Failed to join writer thread " << (writer_thread) <<
           "error code = " << i ;
    PrintMessageDieOnError(msg.str(), WARNING);
  }
#ifdef DEBUG_THREADS
  PrintMessageDieOnError("All thread terminated.", PROGRESS);
#endif

}

int main(int argc, char* argv[]) {
  time_t starttime, processing_starttime,endtime;
  time(&starttime);
  parse_commandline_options(argc, argv);
  if (!quiet) PrintLobSTR();
  unit_name = paired?"pairs":"reads";
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
  // Check that we are using the correct index version
  CheckIndexVersion();
  // Set chrom sizes
  LoadChromSizes();
  // Load referencpe
  LoadReference();

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
  // Set additional alignment params
  opts->max_hits_quit_aln = max_hits_quit_aln;

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

  // run detection/alignment
  PrintMessageDieOnError("Running detection/alignment...", PROGRESS);
  time(&processing_starttime);
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
  time(&endtime);
  run_info.endtime = GetTime();
  DestroyReference();
  OutputRunStatistics();
  OutputRunningTimeInformation(starttime,processing_starttime,endtime,
			       threads, run_info.num_processed_units);
  return 0;
}
