//============================================================================
// Author      : Melissa Gymrek
//============================================================================

#include <boost/algorithm/string.hpp>
#include <err.h>
#include <error.h>
#include <getopt.h>
#include <iostream>
#include <limits.h>
#include <list>
#include <map>
#include <stdlib.h>
#include <string>
#include <unistd.h>

#include "BWAReadAligner.h"
#include "bwtaln.h"
#include "bwase.h"
#include "common.h"
#include "FastaFileReader.h"
#include "FastqFileReader.h"
#include "FFT_nuc_vectors.h"
#include "HammingWindowGenerator.h"
#include "IFileReader.h"
#include "MSReadRecord.h"
#include "MultithreadData.h"
#include "SamFileWriter.h"
#include "STRDetector.h"
#include "runtime_parameters.h"
#include "TabFileWriter.h"
#include "TukeyWindowGenerator.h"

using namespace std;

// list of input files to process from
vector<string> input_files;
vector<string> input_files1;
vector<string> input_files2;

// keep track of # bases so we can calculate coverage
long bases = 0;

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

void show_help(){
  const char* help = "\n\nlobSTR [OPTIONS] {-f <file1[,file2,...]>|--p1 <file1_1[,file2_1,...]> --p2 <file1_2[,file2_1,...]>} --index-prefix <index prefix> -o <output prefix>\n" \
"-f,--files     file or comma-separated list of files containing reads in fasta, fastq, or bam format (default: fasta)\n" \
"--p1             file or comma-separated list of files containing the first end of paired end reads in fasta or fastq (default: fasta)\n" \
"--p2             file or comma-separated list of files containing the second end of paired end reads in fasta or fastq (default: fasta)\n" \
"--gzip           The input files are gzipped (only works for fasta or fastq input)\n" \
"-o,--out       prefix for out put files. will output:\n" \
"                      <prefix>.aligned.tab: tab delimited file of alignments\n" \
"                      <prefix>.aligned.bam: bam file of alignments\n" \
"--index-prefix     prefix for bwa reference (must run lobstr_index.py to create index)\n" \
"                   If the index is downloaded to PATH_TO_INDEX, this argument is PATH_TO_INDEX/lobSTR_\n" \

"\n\nOptions:\n" \
"-h,--help                  display this help screen\n" \
"-v,--verbose               print out useful progress messages\n" \
"-q,--fastq                 reads are in fastq format (default: fasta)\n" \
"--bam                      reads are in bam format (default: fasta)\n" \
"--bampair                  reads are in bam format and are paired-end\n" \
"-p,--threads <INT>         number of threads (default: 1)\n" \
"\n\nAdvanced options - general:\n" \
"--min-read-length <INT>    minimum number of nucleotides for a read to be processed. This should be at least two times fft-window-size (default: 45)\n" \
"--max-read-length <INT>   maximum number of nucleotides for a read to be processed. (default: 1024)\n" \
"\n\nAdvanced options - detection:\n" \
"--fft-window-size    size of fft window (default: 24)\n" \
"--fft-window-step    step size of sliding window (default: 12)\n" \
"--entropy-threshold     threshold score to call a window periodic (defualt: 0.45)\n" \
"--minperiod          minimum period to attempt to detect (default: 2)\n" \
"--maxperiod          maximum period to attempt to detect (default: 6)\n" \
"--minflank <INT>          minimum length of flanking region to try to align (default: 10)\n" \
"--maxflank <INT>          length to trim the ends of flanking regions to if they exceed that length (default: 25)\n" \
"--extend-flank  <INT>     length to extend flanking regions (default: 6)\n" \
"\n\nAdvanced options - alignment:\n" \
"--max-diff-ref       maximum difference in length from the reference sequence to allow for alignment (default: 30) (will take the absolute value)\n" \
"--nw-score           minimum required smith waterman score (maximum is 2*read length)\n" \
"-u                   require length difference to be a multiple of the repeat unit\n" \
    "-m,--mismatch <int>  edit distance allowed during alignment of each flanking region (default: -1)\n" \
    "-g <int>             maximum number of gap opens allowed in each flanking region (default: 1)\n" \
    "-e <int>             maximum number of gap extensions allowed in each flanking region (default: 1)\n" \
    "-r <float>           edit distance allowed during alignment of each flanking region (ignored if -m is set) (default: 0.01)\n" \
"This program takes in raw reads, detects and aligns reads containing microsatellites, and genotypes STR locations.\n\n";
	cout << help;
	exit(1);
}

/*
 * parse the command line options
 */
void parse_commandline_options(int argc,char* argv[]) {
	enum LONG_OPTIONS{
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
	  OPT_SAM,
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
	  OPT_NOTAB,
	  OPT_PARTIALDEBUG,
	  OPT_ORIG_READ,
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
	  {"debug-adjust", 0, 0, OPT_DEBUGADJUST},
	  {"why-not", 0, 0, OPT_WHY_NOT},
	  {"no-tab", 0, 0, OPT_NOTAB},
	  {"partial-debug",0,0,OPT_PARTIALDEBUG},
	  {"orig",0,0,OPT_ORIG_READ},
	  {NULL, no_argument, NULL, 0},
	};
	ch = getopt_long(argc,argv,"hvqp:f:t:g:o:m:s:d:e:g:r:u?",
			 long_options,&option_index);
	while (ch != -1) { 
	  switch(ch){
	  case 'u':
	    unit++;
	    user_defined_arguments += "unit=True;";
	    break;
	  case OPT_PROFILE:
	    profile++;
	    user_defined_arguments += "profile=True;";
	    break;
	  case OPT_PARTIALDEBUG:
	    partial_debug++;
	    break;
	  case OPT_WHY_NOT:
	    why_not_debug++;
	    break;
	  case OPT_ALIGN_DEBUG:
	    align_debug ++;
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
	    break;
	  case 'q':
	  case OPT_FASTQ:
	    input_type = INPUT_FASTQ;
	    fastq++;
	    user_defined_arguments += "fastq=True;";
	    break;
	  case OPT_BAM:
	    input_type = INPUT_BAM;
	    user_defined_arguments += "bam=True;";
	    bam++;
	    break;
	  case OPT_BAMPAIR:
	    paired = true;
	    user_defined_arguments += "bampair=True;";
	    break;
	  case 'p':
	  case OPT_THREADS:
	    threads = atoi(optarg);
	    user_defined_arguments += "threads=";
	    user_defined_arguments += string(optarg);
	    user_defined_arguments += ";";
	    if (threads <=0)
	      errx(1,"Error: invalid number of threads");
	    break;
	  case 'f':
	  case OPT_FILES:
	    input_files_string = optarg;
	    user_defined_arguments += "files=";
	    user_defined_arguments += input_files_string;
	    user_defined_arguments += ";";
	    break;
	  case OPT_PAIR1:
	    input_files_string_p1 = optarg;
	    user_defined_arguments += "files1=";
	    user_defined_arguments += input_files_string_p1;
	    user_defined_arguments += ";";
	    paired = true;
	    break;
	  case OPT_PAIR2:
	    input_files_string_p2 = optarg;
	    user_defined_arguments += "files2=";
	    user_defined_arguments += input_files_string_p2;
	    user_defined_arguments += ";";
	    paired = true;
	    break;
	  case OPT_GZIP:
	    user_defined_arguments += "input_gzipped;";
	    gzip++;
	    break;
	  case 'o':
	  case OPT_OUTPUT:
	    output_prefix = string(optarg);
	    user_defined_arguments += "out=";
	    user_defined_arguments += output_prefix;
	    user_defined_arguments += ";";
	    break;
	  case 'm':
	  case OPT_MISMATCH:
	    allowed_mismatches = atoi(optarg);
	    if (allowed_mismatches < 0)
	      errx(1,"Error: invalid number of mismatches");
	    user_defined_arguments += "m=";
	    user_defined_arguments += string(optarg);
	    user_defined_arguments += ";";
	    break;
	  case 'b':
	  case OPT_SAM:
	    sam++;
	    break;
	  case OPT_FFT_WINDOW_SIZE:
	    fft_window_size = atoi(optarg);
	    user_defined_arguments += "fft-window-size=";
	    user_defined_arguments += string(optarg);
	    user_defined_arguments += ";";
	    break;
	  case OPT_FFT_WINDOW_STEP:
	    fft_window_step = atoi(optarg);
	    user_defined_arguments += "fft-window-step=";
	    user_defined_arguments += string(optarg);
	    user_defined_arguments += ";";
	    break;
	  case OPT_EXTEND:
	    extend = atoi(optarg);
	    if (extend <=0)
	      errx(1, "Error: invalid extension length");
	    user_defined_arguments += "extend=";
	    user_defined_arguments += string(optarg);
	    user_defined_arguments += ";";
	    break;
	  case OPT_MIN_PERIOD:
	    min_period = atoi(optarg);
	    if (min_period <= 0)
	      errx(1, "Error: invalid min period");
	    user_defined_arguments += "minperiod=";
	    user_defined_arguments += string(optarg);
	    user_defined_arguments += ";";
	    break;
	  case OPT_MAX_PERIOD:
	    max_period = atoi(optarg);
	    if (max_period <= 0)
	      errx(1, "Error: invalid max period");
	    user_defined_arguments += "maxperiod=";
	    user_defined_arguments += string(optarg);
	    user_defined_arguments += ";";
	    break;
	  case OPT_MAX_FLANK_LEN:
	    max_flank_len = atoi(optarg);
	    if (max_flank_len <= 0)
	      errx(1, "Error: invalid max flank length");
	    user_defined_arguments += "maxflank=";
	    user_defined_arguments += string(optarg);
	    user_defined_arguments += ";";
	    break;
	  case OPT_MIN_FLANK_LEN:
	    min_flank_len = atoi(optarg);
	    if (min_flank_len <= 0)
	      errx(1, "Error: invalid min flank length");
	    user_defined_arguments += "minflank=";
	    user_defined_arguments += string(optarg);
	    user_defined_arguments += ";";
	    break;
	  case OPT_MAX_DIFF_REF:
	    max_diff_ref = atoi(optarg);
	    if (max_diff_ref <=0 )
	      errx(1, "Error: invalid max diff ref");
	    user_defined_arguments += "max-diff-ref=";
	    user_defined_arguments += string(optarg);
	    user_defined_arguments += ";";
	    break;
	  case OPT_FFTW_DEBUG:
	    fftw_debug = true;
	    break;
	  case OPT_LOBE_DEBUG:
	    lobe_debug = true;
	    break;
	  case OPT_MIN_READ_LENGTH:
	    min_read_length = atoi(optarg);
	    user_defined_arguments += "min-read-length=";
	    user_defined_arguments += string(optarg);
	    user_defined_arguments += ";";
	    break;
	  case OPT_MAX_READ_LENGTH:
	    max_read_length = atoi(optarg);
	    user_defined_arguments += "max-read-length=";
	    user_defined_arguments += string(optarg);
	    user_defined_arguments += ";";
	    break;
	  case OPT_ENTROPY_THRESHOLD:
	    entropy_threshold = atof(optarg);
	    user_defined_arguments += "entropy-threshold=";
	    user_defined_arguments += string(optarg);
	    user_defined_arguments += ";";
	    break;
	  case OPT_DEBUG_ENTROPY:
	    entropy_debug++;
	    break;
	  case OPT_INDEX:
	    index_prefix = string(optarg);
	    user_defined_arguments += "index-prefix=";
	    user_defined_arguments += index_prefix;
	    user_defined_arguments += ";";
	    break;
	  case OPT_GAP_OPEN:
	  case 'g':
	    gap_open = atoi(optarg);
	    user_defined_arguments += "gap-open=";
	    user_defined_arguments += string(optarg);
	    user_defined_arguments += ";";
	    break;
	  case OPT_GAP_EXTEND:
	  case 'e':
	    gap_extend = atoi(optarg);
	    user_defined_arguments += "gap-extend=";
	    user_defined_arguments += string(optarg);
	    user_defined_arguments += ";";
	    break;
	  case OPT_FPR:
	  case 'r':
	    fpr = atof(optarg);
 	    user_defined_arguments += "fpr=";
	    user_defined_arguments += string(optarg);
	    user_defined_arguments += ";";
	    break;
	  case OPT_EXTEND_FLANK:
	    extend_flank = atoi(optarg);
	    user_defined_arguments += "extend-flank=";
	    user_defined_arguments += string(optarg);
	    user_defined_arguments += ";";
	    break;
	  case OPT_SW:
	    min_sw_score = atoi(optarg);
	    user_defined_arguments += "min-sw-score=";
	    user_defined_arguments += string(optarg);
	    user_defined_arguments += ";";
	    break;
	  case OPT_NOTAB:
	    notab++;
	    break;
	  case OPT_ORIG_READ:
	    include_orig_read_start++;
	    user_defined_arguments += "orig-read-start=True;";
	    break;
	  case '?':
	    show_help();
	    exit(1);
	    break;
	  default:
	    show_help();
	    exit(1);
	  }
	  ch = getopt_long(argc,argv,"hvqp:f:t:g:o:m:s:d:e:g:r:u?",
			   long_options,&option_index);

	} 

	// any arguments left over are extra
	if (optind < argc) {
	  cout << "Unnecessary leftover arguments...\n";
	  show_help();
	  exit(1);
	}
	// make sure arguments make sense
	if (fft_window_step > fft_window_size) {
	  errx(1, "fft_window_step must be <=fft_window_size");
	}
	if (min_period > max_period) {
	  errx(1, "min_period must be <= max_period");
	}
	if (min_flank_len > max_flank_len) {
	  errx(1, "min_flank_len must be <=max_flank_len");
	}
	if (min_period < 2 || max_period > 6) {
	  errx(1, "lobSTR can currently only profile STRs of periods 2 through 6.");
	}
	// check that we have the mandatory parameters
	if ( (((!paired || bam) && input_files_string.empty()) || 
	      (paired && !bam && (input_files_string_p1.empty() || input_files_string_p2.empty())))||  
	      output_prefix.empty() || index_prefix.empty()) {
	  errx(1, "Required arguments are mising");
	} 
	if (gzip && bam ) {
	  errx(1, "Gzip option not compatible with bam input");
	}
}

void LoadReference(const std::string& repseq) {
  // Load BWT index
  char* prefix = (char*)calloc(strlen(index_prefix.c_str()) +
			       strlen(repseq.c_str()) + 20, 1);
  strcpy(prefix, index_prefix.c_str()); strcat(prefix, repseq.c_str());
  strcat(prefix, ".fa");
  //  cout << "loading ref " << repseq << endl;
  bwt_t *bwt[2];
  char *str = (char*)calloc(strlen(prefix) + 10, 1);
  strcpy(str, prefix); strcat(str, ".bwt");
  bwt[0] = bwt_restore_bwt(str);
  strcpy(str, prefix); strcat(str, ".rbwt");
  bwt[1] = bwt_restore_bwt(str);
  strcpy(str, prefix); strcat(str, ".sa");
  bwt_restore_sa(str, bwt[0]);
  strcpy(str, prefix); strcat(str, ".rsa");
  bwt_restore_sa(str, bwt[1]);

  BWT bwt_ref;
  bwt_ref.bwt[0] = bwt[0];
  bwt_ref.bwt[1] = bwt[1];
  bwt_references.insert(pair<string, BWT>(repseq,bwt_ref));

  // Load BNT annotations
  BNT bnt;
  bntseq_t *bns;
  bns = bns_restore(prefix);
  bnt.bns = bns;
  bnt_annotations.insert(pair<string, BNT>(repseq,bnt));
  free(str);
  free(prefix);
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
  TabFileWriter pWriter(output_prefix + ".aligned.tab");
  SamFileWriter samWriter(output_prefix + ".aligned.bam", chrom_sizes);
  STRDetector *pDetector = new STRDetector();
  BWAReadAligner *pAligner = new BWAReadAligner(&bwt_references, &bnt_annotations, &ref_sequences, opts);
  std::string file1;
  std::string file2;
  for (size_t i = 0; i < input_files1.size(); i++) {
    file1 = input_files1.at(i);
    if (paired && !bam) {
      file2 = input_files2.at(i);
      if (!(fexists(file1.c_str()) && fexists(file2.c_str()))) {
	cerr << "Warning: file " << file1 << " or " << file2 << " does not exist" << endl;
	continue;
      }
    } else {
      if (!fexists(file1.c_str())) {
	cerr << "Warning: file " << file1 << " does not exist" << endl;
	continue;
      }
    }
    if (paired) {
      cout << "processing files " << file1 << " and " << file2 << "...\n";
    } else {
      cout << "processing file " <<  file1 << " ...\n";
    }
    IFileReader* pReader = create_file_reader(file1, file2);
    int aligned = false;
    int num_reads_processed = 0;
    std::string repseq = "";
    while (pReader->GetNextRecord(&read_pair)) {
      aligned = false;
      num_reads_processed += 1;
      read_pair.read_count = num_reads_processed;
      // reset fields
      read_pair.reads.at(0).msRepeat = "";
      if (paired) read_pair.reads.at(1).msRepeat = "";
      // Check read length
      if (!(read_pair.reads.at(0).nucleotides.length() >= min_read_length) &&
	  (read_pair.reads.at(0).nucleotides.length() <= max_read_length)) {
	continue;
      }
      if (paired) {
	if (!(read_pair.reads.at(1).nucleotides.length() >= min_read_length) &&
	    (read_pair.reads.at(1).nucleotides.length() <= max_read_length)) {
	  continue;
	}
      }
      bases += read_pair.reads.at(0).nucleotides.length();
      if (paired) bases += read_pair.reads.at(1).nucleotides.length();

      // STEP 1: Sensing
      if (!pDetector->ProcessReadPair(&read_pair)) {
	continue;
      }

      // STEP 2: Alignment
      if (pAligner->ProcessReadPair(&read_pair)) {
	aligned = true;
      } else {
	read_pair.read1_passed_detection = false;
	read_pair.read2_passed_detection = false;
	// Try second best period for each read
	if (read_pair.reads.at(0).ms_repeat_next_best_period != 0) {
	  read_pair.reads.at(0).ms_repeat_best_period =
	    read_pair.reads.at(0).ms_repeat_next_best_period;
	  if (getMSSeq(read_pair.reads.at(0).detected_ms_region_nuc,
		       read_pair.reads.at(0).ms_repeat_best_period, &repseq)) {
	    read_pair.reads.at(0).repseq = repseq;
	    read_pair.read1_passed_detection = true;
	  }
	}
	if (paired) {
	  if (read_pair.reads.at(1).ms_repeat_next_best_period != 0) {
	    read_pair.reads.at(1).ms_repeat_best_period =
	      read_pair.reads.at(1).ms_repeat_next_best_period;
	    if (getMSSeq(read_pair.reads.at(1).detected_ms_region_nuc,
			 read_pair.reads.at(1).ms_repeat_best_period, &repseq)) {
	      read_pair.reads.at(1).repseq = repseq;
	      read_pair.read2_passed_detection = true;
	    }
	  }
	}
	if (read_pair.read1_passed_detection ||
	    read_pair.read2_passed_detection) {
	  if (pAligner->ProcessReadPair(&read_pair)) {
	    aligned = true;
	  }
	}
      }
      if (aligned) {
	pWriter.WriteRecord(read_pair);
	if (sam) samWriter.WriteRecord(read_pair);
      }
    }
    delete pReader;
    cout << "Processed " << num_reads_processed << " reads" << endl;
  }
  delete pDetector;
  delete pAligner;
}


void* satellite_process_consumer_thread(void *arg) {
  MultithreadData *pMT_DATA = (MultithreadData*)(arg);
  STRDetector *pDetector = new STRDetector();
  BWAReadAligner *pAligner = new BWAReadAligner(&bwt_references, &bnt_annotations, &ref_sequences, opts);
  int aligned = false;
  while (1) {
    aligned = false;
    std::string repseq;
    ReadPair* pReadRecord = pMT_DATA->get_new_input();
    if (!(pReadRecord->reads.at(0).nucleotides.length() >= min_read_length) &&
	(pReadRecord->reads.at(0).nucleotides.length() <= min_read_length)) {
      continue;
    }
    if (paired) {
      if (!(pReadRecord->reads.at(1).nucleotides.length() >= min_read_length) &&
	  (pReadRecord->reads.at(1).nucleotides.length() <= min_read_length)) {
	continue;
      }
    }
    bases += pReadRecord->reads.at(0).nucleotides.length();
    if (paired) pReadRecord->reads.at(1).nucleotides.length();
    
    // STEP 1: Sensing
    if (!pDetector->ProcessReadPair(pReadRecord)) {
      continue;
    }

    // STEP 2: Alignment
    if (pAligner->ProcessReadPair(pReadRecord)) {
      aligned = true;
    } else {
      pReadRecord->read1_passed_detection = false;
      pReadRecord->read2_passed_detection = false;
      // Try second best period for each read
      if (pReadRecord->reads.at(0).ms_repeat_next_best_period != 0) {
	pReadRecord->reads.at(0).ms_repeat_best_period =
	  pReadRecord->reads.at(0).ms_repeat_next_best_period;
	if (getMSSeq(pReadRecord->reads.at(0).detected_ms_region_nuc,
		     pReadRecord->reads.at(0).ms_repeat_best_period, &repseq)) {
	  pReadRecord->reads.at(0).repseq = repseq;
	  pReadRecord->read1_passed_detection = true;
	}
      }
      if (paired) {
	if (pReadRecord->reads.at(1).ms_repeat_next_best_period != 0) {
	  pReadRecord->reads.at(1).ms_repeat_best_period =
	    pReadRecord->reads.at(1).ms_repeat_next_best_period;
	  if (getMSSeq(pReadRecord->reads.at(1).detected_ms_region_nuc,
		       pReadRecord->reads.at(1).ms_repeat_best_period, &repseq)) {
	    pReadRecord->reads.at(1).repseq = repseq;
	    pReadRecord->read2_passed_detection = true;
	  }
	}
      }
      if (pReadRecord->read1_passed_detection ||
	  pReadRecord->read2_passed_detection) {
	if (pAligner->ProcessReadPair(pReadRecord)) {
	  aligned = true;
	}
      }
    }
    if (aligned) {
      pMT_DATA->post_new_output_read(pReadRecord);
    } else {
      //Don't pass this read on to the output thread,
      //discard it here
      delete pReadRecord;
      pMT_DATA->increment_output_counter();
    }
  }
  pthread_exit((void*) arg);
}

void* output_writer_thread(void *arg) {
  MultithreadData *pMT_DATA = (MultithreadData*)(arg);
  TabFileWriter *pWriter = new TabFileWriter(output_prefix + ".aligned.tab");
  SamFileWriter samWriter(output_prefix + ".aligned.bam", chrom_sizes);
  while (1) {
    ReadPair *pReadRecord = pMT_DATA->get_new_output();
    pWriter->WriteRecord(*pReadRecord);
    samWriter.WriteRecord(*pReadRecord);
    delete pReadRecord;
    pMT_DATA->increment_output_counter();
  }
  delete pWriter;
  delete pMT_DATA;
  pthread_exit((void*) arg);
}

void multi_thread_process_loop(vector<string> files1, vector<string> files2) {
  MultithreadData mtdata(threads);
  list<pthread_t> satellite_threads;
  pthread_t writer_thread;
  if (files1.size() == 0) return;
  for (size_t i=0; i<threads; ++i) {
    pthread_t id;
    if (pthread_create(&id, NULL, satellite_process_consumer_thread, (void*)&mtdata))
      err(1,"Failed to create thread");
    satellite_threads.push_back(id);
  }
  
  if (pthread_create(&writer_thread, NULL, output_writer_thread, (void*)&mtdata))
    err(1,"failed to create output writer thread");
  
  size_t counter = 1 ;
  std::string file1;
  std::string file2;
  for (size_t i = 0; i < files1.size(); i++) {
    file1 = input_files1.at(i);
    if (paired && !bam) {
      file2 = input_files2.at(i);
      if (!(fexists(file1.c_str()) && fexists(file2.c_str()))) {
	cerr << "Warning: file " << file1 << " or " << file2 << " does not exist" << endl;
	continue;
      }
    } else {
      if (!fexists(file1.c_str())) {
	cerr << "Warning: file " << file1 << " does not exist" << endl;
	continue;
      }
    }
    if (paired) {
      cout << "processing files " << file1 << " and " << file2 << "...\n";
    } else {
      cout << "processing file " <<  file1 << " ...\n";
    }
    IFileReader *pReader = create_file_reader(file1, file2);
    do {
      ReadPair *pRecord = new ReadPair;
      pRecord->read_count = counter;
      if (!pReader->GetNextRecord(pRecord))
	break; //no more reads
	counter++;
	mtdata.increment_input_counter();	
	mtdata.post_new_input_read(pRecord);
	pRecord = NULL ; //the consumers will take it from here, and free it
      } while(1);
      delete pReader;
  }
  while ( 1 ) {
    sleep(1); //OMG, the horror...
    if ( mtdata.input_output_counters_equal())
      break;
    sleep(10);
    break;
  }
}

int main(int argc,char* argv[]) {
  parse_commandline_options(argc,argv);
  
  if (my_verbose) cout << "Initializing..." << endl;
  // open file with all names
  TextFileReader tReader(index_prefix+"strdict.txt");
  string line;
  while(tReader.GetNextLine(&line)) {
    // set chrom sizes
    if (line.substr(0,3) == "REF") {
      vector<string> items;
      split(line,'\t', items);
      chrom_sizes[items[1]] = atoi(items[2].c_str());
    } else {
      // make sure repeat is valid
      if (count(line, line.at(0)) != line.length()) {
	LoadReference(line);
      }
    }
  }

  // Load fasta reference with all STR sequences
  FastaFileReader faReader(index_prefix+"ref.fa");
  MSReadRecord ref_record;
  while(faReader.GetNextRead(&ref_record)) {
    REFSEQ refseq;
    vector<string> items;
    string refstring = ref_record.ID;;
    split(refstring, '_', items);
    if (items.size() == 7) { // else must be _random_, _cox_hap1, etc.
      refseq.sequence = ref_record.orig_nucleotides;
      refseq.start = atoi(items.at(2).c_str());
      int refid = atoi(items.at(0).c_str());
      ref_sequences.insert(pair<int,REFSEQ>(refid, refseq));
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
  opts->max_gapo=gap_open;
  opts->max_gape=gap_extend;
  // take the first INT subsequence as seed
  opts->seed_len = 5;
  // all hits with no more than maxdiff found

  // get the input files
  if (paired) {					
    boost::split(input_files1, input_files_string_p1, boost::is_any_of(","));
    boost::split(input_files2, input_files_string_p2, boost::is_any_of(","));
    if (!input_files1.size() == input_files2.size()) {
      errx(1, "Error: different number of files for each pair");
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
  if (my_verbose) {cerr << "Running detection/alignment..." << endl;}
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
  return 0;
}
